using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using Random = UnityEngine.Random;

public class SPH : MonoBehaviour
{
    const int NumThreads = 8;

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
    public int gridResolution = 64;

    #endregion

    #region Private

    // Particle
    private float particleRadius;
    private float effectiveRadius;

    // Wall Weight
    private RenderTexture wallWeightTexture;

    // Distance
    private int gridResolutionX, gridResolutionY, gridResolutionZ;
    private RenderTexture distanceTexture;
    private float cellSize;
    private MeshFilter meshFilter;

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

        cellSize = scale.y / gridResolution;

        gridResolutionY = gridResolution;
        gridResolutionX = Mathf.CeilToInt(scale.x / cellSize);
        gridResolutionZ = Mathf.CeilToInt(scale.z / cellSize);

        transform.localScale = new(gridResolutionX * cellSize, scale.y, gridResolutionZ * cellSize);

        distanceTexture = CreateRenderTexture3D(gridResolutionX, gridResolutionY, gridResolutionZ, RenderTextureFormat.RFloat);

        camera = Camera.main;
        cameraOrbit = camera.GetComponent<CameraOrbit>();
        InitializeBoxes();

        // Adjust particle resolution for better distribution
        var particlesWidth = Mathf.NextPowerOfTwo(Mathf.CeilToInt(Mathf.Sqrt(desiredParticleCount)));
        var particlesHeight = Mathf.CeilToInt((float)desiredParticleCount / particlesWidth);

        var particleCount = particlesWidth * particlesHeight;

        // Adjust particle radius based on grid resolution
        Vector3[] particlesPositions = CreateParticlePositions(particleCount);

        UpdateDistanceTexture(meshFilter);
        InitializeWallWeights();
    }

    private void Update()
    {
        UpdateMouse(out Vector3 worldSpaceMouseRay, out Vector3 worldMouseVelocity);

        if (renderParticles)
            Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _particleArgsBuffer);
    }

    private void OnDestroy()
    {
        _particleMeshPropertiesBuffer?.Release();
        _particleArgsBuffer?.Release();
        distanceTexture.Release();
    }

    #endregion

    #region Initializations

    private void InitShaders()
    {
        distanceShader = Resources.Load<ComputeShader>("Distance");
        clearShader = Resources.Load<ComputeShader>("Clear");
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

        distanceShader.SetTexture(0, ShaderIDs.Distance, distanceTexture);
        distanceShader.SetBuffer(0, ShaderIDs.Vertices, vertices);
        distanceShader.SetBuffer(0, ShaderIDs.Triangles, triangles);

        int threadGroupsX = Mathf.CeilToInt((float)gridResolutionX / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)gridResolutionY / NumThreads);
        int threadGroupsZ = Mathf.CeilToInt((float)gridResolutionZ / NumThreads);

        distanceShader.Dispatch(0, threadGroupsX, threadGroupsY, threadGroupsZ);

        triangles.Release();
        vertices.Release();
    }

    #endregion

    #region Texture Helpers

    private void Swap(ref RenderTexture a, ref RenderTexture b)
    {
        (b, a) = (a, b);
    }

    private static RenderTexture CreateRenderTexture2D(int width, int height, RenderTextureFormat format)
    {
        var rt = new RenderTexture(width, height, 0, format)
        {
            enableRandomWrite = true,
            filterMode = FilterMode.Bilinear,
            wrapMode = TextureWrapMode.Clamp
        };
        rt.Create();
        return rt;
    }

    private static RenderTexture CreateRenderTexture3D(int width, int height, int depth, RenderTextureFormat format)
    {
        var rt = new RenderTexture(width, height, 0, format)
        {
            dimension = TextureDimension.Tex3D,
            volumeDepth = depth,
            enableRandomWrite = true,
            filterMode = FilterMode.Bilinear,
            wrapMode = TextureWrapMode.Clamp
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
