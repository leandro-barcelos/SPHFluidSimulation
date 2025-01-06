using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using Random = UnityEngine.Random;

public class SPH : MonoBehaviour
{
    const int NumThreads = 8;
    const int MaxParticlesPerVoxel = 32;

    #region Auxiliary Structures
    private struct MeshProperties
    {
        public Matrix4x4 Mat;
        public Vector4 Color;

        public static int Size()
        {
            return
                sizeof(float) * 4 * 4 + // matrix;
                sizeof(float) * 4;      // color;
        }
    }
    #endregion

    #region Public

    [Header("Parameters")]
    [Range(0f, 1f / 60f)] public float timeStep = 1f / 60f;
    public int desiredParticleCount = 1000;
    [Range(0.01f, 1f)] public float spawnAreaFillRate;
    [Range(3, 1000)] public int wallWeightNumSamples = 100;

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial;
    public bool renderParticles;

    [Header("Distance")]
    public int fastSweepingPasses = 3;

    #endregion

    #region Private

    // Particle
    private float particleRadius;
    private float effectiveRadius;
    private RenderTexture particlePositionTexture;
    private int particleWidth;
    private int particleHeight;

    // Wall Weight
    private RenderTexture wallWeightTexture;

    // Distance
    private int gridResolutionX, gridResolutionY, gridResolutionZ;
    private RenderTexture distanceTexture;
    private float cellSize;
    private MeshFilter meshFilter;

    // Bucket
    private ComputeBuffer bucketBuffer;
    private ComputeShader bucketShader;

    // Mouse
    private Vector3 lastMousePlane = Vector3.zero;

    // Camera
    private CameraOrbit cameraOrbit;
    private new Camera camera;

    // Boxes
    private List<Bounds> boxes;

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _particleMeshPropertiesBuffer, _particleArgsBuffer;
    private Bounds _bounds;

    // Shaders
    private ComputeShader distanceShader, clearShader;

    #endregion

    #region Unity Functions

    private void Start()
    {
        InitShaders();

        meshFilter = gameObject.GetComponent<MeshFilter>();

        var scale = transform.localScale;

        camera = Camera.main;
        cameraOrbit = camera.GetComponent<CameraOrbit>();
        InitializeBoxes();

        // Adjust particle resolution for better distribution
        particleWidth = Mathf.NextPowerOfTwo(Mathf.CeilToInt(Mathf.Sqrt(desiredParticleCount)));
        particleHeight = Mathf.CeilToInt((float)desiredParticleCount / particleWidth);

        var particleCount = particleWidth * particleHeight;

        InitializeParticleTextures(particleCount);

        cellSize = 8 * particleRadius;

        gridResolutionX = Mathf.CeilToInt(scale.x / cellSize);
        gridResolutionY = Mathf.CeilToInt(scale.y / cellSize);
        gridResolutionZ = Mathf.CeilToInt(scale.z / cellSize);

        transform.localScale = new(gridResolutionX * cellSize, scale.y, gridResolutionZ * cellSize);

        distanceTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RFloat, FilterMode.Bilinear);

        // Initialize bucket buffer
        int totalBucketSize = gridResolutionX * gridResolutionY * gridResolutionZ * MaxParticlesPerVoxel;
        bucketBuffer = new ComputeBuffer(totalBucketSize, sizeof(uint));

        UpdateDistanceTexture(meshFilter);
        InitializeWallWeights();
    }

    private void Update()
    {
        UpdateMouse(out Vector3 worldSpaceMouseRay, out Vector3 worldMouseVelocity);
        BucketGeneration();

        if (renderParticles)
            Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _particleArgsBuffer);
    }

    private void OnDestroy()
    {
        _particleMeshPropertiesBuffer?.Release();
        _particleArgsBuffer?.Release();
        distanceTexture.Release();
        bucketBuffer?.Release();
        particlePositionTexture.Release();
    }

    #endregion

    #region Initializations

    private void InitShaders()
    {
        distanceShader = Resources.Load<ComputeShader>("Distance");
        clearShader = Resources.Load<ComputeShader>("Clear");
        bucketShader = Resources.Load<ComputeShader>("Bucket");
    private void InitializeParticleTextures(int particleCount)
    {
        // Create particle position texture
        particlePositionTexture = CreateRenderTexture2D(particleWidth, particleHeight, RenderTextureFormat.ARGBFloat);
        particleDensityTexture = CreateRenderTexture2D(particleWidth, particleHeight, RenderTextureFormat.RFloat);

        // Adjust particle radius based on grid resolution
        Vector3[] particlesPositions = CreateParticlePositions(particleCount);
        Texture2D particlePositionTexture2D = new(particleWidth, particleHeight, TextureFormat.RGBAFloat, false);
        Color[] particleColors = new Color[particleCount];

        for (int i = 0; i < particleCount; i++)
        {
            Vector3 pos = particlesPositions[i];
            particleColors[i] = new Color(pos.x, pos.y, pos.z, 1.0f);
        }

        particlePositionTexture2D.SetPixels(particleColors);
        particlePositionTexture2D.Apply();

        Graphics.Blit(particlePositionTexture2D, particlePositionTexture);
        Destroy(particlePositionTexture2D);
    }

    private Vector3[] CreateParticlePositions(int particleCount)
    {
        if (particleCount == 0)
            return new Vector3[0];

        _particleMesh = OctahedronSphereCreator.Create(1, 1f);
        _bounds = new Bounds(transform.position, Vector3.one * (occlusionRange + 1));

        uint[] args = { 0, 0, 0, 0, 0 };
        args[0] = _particleMesh.GetIndexCount(0);
        args[1] = (uint)particleCount;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _particleArgsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _particleArgsBuffer.SetData(args);

        _particleMeshPropertiesBuffer = new ComputeBuffer(particleCount, MeshProperties.Size());

        Vector3[] particlesPositions = new Vector3[particleCount];
        MeshProperties[] properties = new MeshProperties[particleCount];
        var totalVolume = 0f;

        Quaternion rotation = Quaternion.identity;

        // Calculate total volume of spawn boxes
        foreach (var box in boxes)
        {
            var volume = box.size.x * box.size.y * box.size.z;
            totalVolume += volume;
        }

        var particleVol = totalVolume / desiredParticleCount;
        particleRadius = Mathf.Pow(3 * particleVol / (4 * Mathf.PI), 1f / 3f) * spawnAreaFillRate;
        effectiveRadius = 4 * particleRadius;

        Vector3 particleScale = new(particleRadius * 2, particleRadius * 2, particleRadius * 2);

        // Distribute particles among boxes
        var particlesCreatedSoFar = 0;
        for (var i = 0; i < boxes.Count; i++)
        {
            var box = boxes[i];
            var volume = box.size.x * box.size.y * box.size.z;

            int particlesInBox;
            if (i < boxes.Count - 1)
            {
                particlesInBox = Mathf.FloorToInt(particleCount * volume / totalVolume);
            }
            else
            {
                particlesInBox = particleCount - particlesCreatedSoFar;
            }

            // Create particles with jittered positions
            for (var j = 0; j < particlesInBox; j++)
            {
                var position = new Vector3(
                    Random.Range(box.min.x, box.max.x),
                    Random.Range(box.min.y, box.max.y),
                    Random.Range(box.min.z, box.max.z)
                );

                MeshProperties props = new()
                {
                    Mat = Matrix4x4.TRS(position, rotation, particleScale),
                    Color = Color.blue
                };

                particlesPositions[particlesCreatedSoFar + j] = position;
                properties[particlesCreatedSoFar + j] = props;
            }

            particlesCreatedSoFar += particlesInBox;
        }

        _particleMeshPropertiesBuffer.SetData(properties);
        particleMaterial.SetBuffer("_Properties", _particleMeshPropertiesBuffer);

        return particlesPositions;
    }

    private void InitializeBoxes()
    {
        boxes = new List<Bounds>();
        foreach (var meshFilter in FindObjectsByType<MeshFilter>(FindObjectsSortMode.None))
        {
            if (meshFilter.gameObject.CompareTag("Spawn"))
            {
                var center = meshFilter.transform.position;
                var size = meshFilter.transform.localScale;
                var box = new Bounds(center, size);
                boxes.Add(box);
                meshFilter.gameObject.SetActive(false);
            }
        }
    }

    private void InitializeWallWeights()
    {
        var wallParticles = InitializeWallParticles();

        Color[] weights = new Color[wallWeightNumSamples];

        var weightConstant = 315f / (64 * Mathf.PI * Mathf.Pow(effectiveRadius, 9));
        var sqrEffectiveRadius = effectiveRadius * effectiveRadius;

        var step = effectiveRadius / wallWeightNumSamples;
        var position = new Vector3(effectiveRadius, effectiveRadius, effectiveRadius);

        for (var i = 0; i < wallWeightNumSamples; i++)
        {
            foreach (var particle in wallParticles)
            {
                var distance = Vector3.Distance(position, particle);

                if (distance <= effectiveRadius)
                {
                    weights[i] += new Color(weightConstant * Mathf.Pow(sqrEffectiveRadius - (distance * distance), 3), 0f, 0f);
                }
            }

            position.y += step;
        }

        wallWeightTexture = CreateRenderTexture2D(wallWeightNumSamples, 1, RenderTextureFormat.RFloat);
        Texture2D wallWeightTexture2D = CreateTexture2D(wallWeightNumSamples, 1, weights, TextureFormat.RFloat);

        Graphics.Blit(wallWeightTexture2D, wallWeightTexture);
        Destroy(wallWeightTexture2D);
    }

    private List<Vector3> InitializeWallParticles()
    {
        Bounds wallBounds = new(
            new Vector3(effectiveRadius, effectiveRadius / 2f, effectiveRadius),
            new Vector3(2 * effectiveRadius, effectiveRadius, 2 * effectiveRadius) + new Vector3(particleRadius / 2, particleRadius / 2, particleRadius / 2)
        );

        var numParticles = Mathf.FloorToInt(wallBounds.size.x / (particleRadius * 2));

        var spaceBetween = (wallBounds.size.x - numParticles * (particleRadius * 2)) / (numParticles - 1);

        List<Vector3> wallParticles = new List<Vector3>();

        for (int i = 0; i < numParticles; i++)
        {
            for (int j = 0; j < numParticles / 2; j++)
            {
                for (int k = 0; k < numParticles; k++)
                {
                    Vector3 position = new(
                        wallBounds.min.x + i * (particleRadius * 2 + spaceBetween),
                        wallBounds.max.y - j * (particleRadius * 2 + spaceBetween),
                        wallBounds.min.z + k * (particleRadius * 2 + spaceBetween)
                    );

                    if (wallBounds.Contains(position))
                    {
                        wallParticles.Add(position);
                    }
                }
            }
        }

        return wallParticles;
    }

    #endregion

    #region Update

    private void UpdateMouse(out Vector3 worldSpaceMouseRay, out Vector3 worldMouseVelocity)
    {
        worldSpaceMouseRay = Vector3.zero;
        worldMouseVelocity = Vector3.zero;

        if (Input.GetMouseButton(2))
        {
            Ray ray = camera.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(ray, out _))
            {
                Vector3 currentMousePlane = ray.GetPoint(cameraOrbit.distance);
                worldMouseVelocity = (currentMousePlane - lastMousePlane) / Time.deltaTime * 0.6f;
                lastMousePlane = currentMousePlane;
                worldSpaceMouseRay = ray.direction;
            }
        }
    }

    private void UpdateDistanceTexture(MeshFilter meshFilter)
    {
        ClearTexture(distanceTexture);

        int triangleCount = meshFilter.mesh.triangles.Length;
        ComputeBuffer triangles = new(triangleCount, sizeof(int));
        triangles.SetData(meshFilter.mesh.triangles);
        ComputeBuffer vertices = new(meshFilter.mesh.vertexCount, sizeof(float) * 3);
        vertices.SetData(meshFilter.mesh.vertices);

        distanceShader.SetMatrix(ShaderIDs.MeshTransform, meshFilter.transform.localToWorldMatrix);
        distanceShader.SetVector(ShaderIDs.SimOrigin, transform.position - transform.localScale / 2);
        distanceShader.SetVector(ShaderIDs.GridResolution, new(gridResolutionX, gridResolutionY, gridResolutionZ));
        distanceShader.SetFloat(ShaderIDs.CellSize, cellSize);
        distanceShader.SetInt(ShaderIDs.TriangleCount, triangleCount);

        distanceShader.SetTexture(0, ShaderIDs.DistanceTexture, distanceTexture);
        distanceShader.SetBuffer(0, ShaderIDs.Vertices, vertices);
        distanceShader.SetBuffer(0, ShaderIDs.Triangles, triangles);

        int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);

        distanceShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);

        triangles.Release();
        vertices.Release();
    }

    private void BucketGeneration()
    {
        // Clear bucket buffer with particle count as empty marker
        int totalBucketSize = gridResolutionX * gridResolutionY * gridResolutionZ * MaxParticlesPerVoxel;
        uint[] clearData = new uint[totalBucketSize];
        for (int i = 0; i < totalBucketSize; i++)
            clearData[i] = (uint)(particleWidth * particleHeight);
        bucketBuffer.SetData(clearData);

        // Set shader parameters
        bucketShader.SetBuffer(0, ShaderIDs.Bucket, bucketBuffer);
        bucketShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTexture);
        bucketShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleWidth, particleHeight));
        bucketShader.SetVector(ShaderIDs.BucketResolution, new(gridResolutionX, gridResolutionY, gridResolutionZ));
        bucketShader.SetVector(ShaderIDs.SimOrigin, transform.position - transform.localScale / 2);
        bucketShader.SetVector(ShaderIDs.SimScale, transform.localScale);

        // Calculate dispatch size for 2D thread groups
        int threadGroupsX = Mathf.CeilToInt((float)particleWidth / 32);
        int threadGroupsY = Mathf.CeilToInt((float)particleHeight / 32);
        bucketShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
    }

    #endregion

    #region Texture Helpers

    private void Swap(ref RenderTexture a, ref RenderTexture b)
    {
        (b, a) = (a, b);
    }

    private static RenderTexture CreateRenderTexture2D(int width, int height, RenderTextureFormat format, FilterMode filterMode = FilterMode.Point, TextureWrapMode wrapMode = TextureWrapMode.Clamp)
    {
        var rt = new RenderTexture(width, height, 0, format)
        {
            enableRandomWrite = true,
            filterMode = filterMode,
            wrapMode = wrapMode
        };
        rt.Create();
        return rt;
    }

    private static RenderTexture CreateRenderTexture3D(int width, int height, int depth, RenderTextureFormat format, FilterMode filterMode = FilterMode.Point, TextureWrapMode wrapMode = TextureWrapMode.Clamp)
    {
        var rt = new RenderTexture(width, height, 0, format)
        {
            dimension = TextureDimension.Tex3D,
            volumeDepth = depth,
            enableRandomWrite = true,
            filterMode = filterMode,
            wrapMode = wrapMode
        };
        rt.Create();
        return rt;
    }

    private static Texture2D CreateTexture2D(int width, int height, Color[] colors, TextureFormat format)
    {
        var texture = new Texture2D(width, height, format, false);
        texture.SetPixels(colors);
        texture.Apply();
        return texture;
    }

    private void ClearTexture(RenderTexture texture)
    {
        clearShader.SetTexture(0, ShaderIDs.Texture, texture);
        int threadGroupsX = Mathf.CeilToInt((float)texture.width / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)texture.height / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)texture.volumeDepth / NumThreads);
        clearShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);
    }

    #endregion
}
