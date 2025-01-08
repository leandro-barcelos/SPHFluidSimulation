using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using Random = UnityEngine.Random;

public class SPH : MonoBehaviour
{
    #region Constants

    const int NumThreads = 8;
    const int MaxParticlesPerVoxel = 32;
    const int Read = 0;
    const int Write = 1;

    #endregion

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

    [Header("Particle Initialization")]
    [Range(1024, 4194304)] public int desiredParticleCount = 1024;
    [Range(1, 256)] public int bucketResolution = 256;
    [Range(0.01f, 1f)] public float initialFluidHeight;

    [Header("Parameters")]
    [Range(0f, 1f / 60f)] public float timeStep = 1f / 60f;
    [Range(0.001f, 0.01f)] public float viscosity = 0.01f;
    [Range(1f, 2f)] public float restDensity = 1.5f;
    [Range(100f, 500f)] public float gasConstant = 150.0f;
    [Range(1000f, 10000f)] public float stiffnessCoeff = 5000.0f;
    [Range(1f, 50f)] public float dampingCoeff = 10.0f;

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial;
    public bool renderParticles;

    #endregion

    #region Private

    // Particle
    private float effectiveRadius;
    private RenderTexture[] particlePositionTextures;
    private RenderTexture[] particleVelocityTextures;
    private RenderTexture particleDensityTexture;
    private int particleTextureResolution;
    private float particleMass;

    // Bucket
    private ComputeBuffer bucketBuffer;

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _particleMeshPropertiesBuffer, _particleArgsBuffer;
    private Bounds _bounds;

    // Shaders
    private ComputeShader bucketShader, clearShader, densityShader, velPosShader;

    #endregion

    #region Unity Functions

    private void Start()
    {
        desiredParticleCount = Mathf.NextPowerOfTwo(desiredParticleCount);
        particleTextureResolution = (int)Mathf.Sqrt(desiredParticleCount);

        InitShaders();

        CreateParticleTextures();

        InitializeParticles();

        InitializeBucketBuffer();

        InitializeCameraOrbit();
    }

    private void Update()
    {
        BucketGeneration();
        DensityCalculation();
        UpdateVelocityAndPosition();

        if (renderParticles)
            Graphics.DrawMeshInstancedIndirect(_particleMesh, 0, particleMaterial, _bounds, _particleArgsBuffer);
    }

    private void OnDestroy()
    {
        _particleMeshPropertiesBuffer?.Release();
        _particleArgsBuffer?.Release();
        bucketBuffer?.Release();
        particlePositionTextures[Read].Release();
        particleVelocityTextures[Read].Release();
        particlePositionTextures[Write].Release();
        particleVelocityTextures[Write].Release();
        particleDensityTexture?.Release();
    }

    #endregion

    #region Initializations

    private void InitShaders()
    {
        clearShader = Resources.Load<ComputeShader>("Clear");
        bucketShader = Resources.Load<ComputeShader>("Bucket");
        densityShader = Resources.Load<ComputeShader>("Density");
        velPosShader = Resources.Load<ComputeShader>("VelPos");
    }

    private void CreateParticleTextures()
    {
        // Create particle position textures
        particlePositionTextures = new RenderTexture[2];

        particlePositionTextures[Read] = CreateRenderTexture2D(particleTextureResolution, particleTextureResolution, RenderTextureFormat.ARGBFloat);
        particlePositionTextures[Write] = CreateRenderTexture2D(particleTextureResolution, particleTextureResolution, RenderTextureFormat.ARGBFloat);

        // Create particle velocity textures
        particleVelocityTextures = new RenderTexture[2];

        particleVelocityTextures[Read] = CreateRenderTexture2D(particleTextureResolution, particleTextureResolution, RenderTextureFormat.ARGBFloat);
        particleVelocityTextures[Write] = CreateRenderTexture2D(particleTextureResolution, particleTextureResolution, RenderTextureFormat.ARGBFloat);

        // Create density texture
        particleDensityTexture = CreateRenderTexture2D(particleTextureResolution, particleTextureResolution, RenderTextureFormat.RFloat);

    }

    private void InitializeParticles()
    {
        effectiveRadius = 1f / (bucketResolution - 1);
        particleMass = damFillRate / particleNumber;

        // Initialize render properties
        _particleMesh = OctahedronSphereCreator.Create(1, 1f);
        _bounds = new Bounds(transform.position, Vector3.one * (occlusionRange + 1));

        uint[] args = { 0, 0, 0, 0, 0 };
        args[0] = _particleMesh.GetIndexCount(0);
        args[1] = (uint)desiredParticleCount;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _particleArgsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _particleArgsBuffer.SetData(args);

        _particleMeshPropertiesBuffer = new ComputeBuffer(desiredParticleCount, MeshProperties.Size());

        Vector3[] particlesPositions = CreateParticlePositions();

        particleMaterial.SetBuffer(ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        // Set texture to positions
        Texture2D particlePositionTexture2D = new(particleTextureResolution, particleTextureResolution, TextureFormat.RGBAFloat, false);
        Color[] particleColors = new Color[desiredParticleCount];

        for (int i = 0; i < desiredParticleCount; i++)
        {
            Vector3 pos = particlesPositions[i];
            particleColors[i] = new Color(pos.x, pos.y, pos.z, 1.0f);
        }

        particlePositionTexture2D.SetPixels(particleColors);
        particlePositionTexture2D.Apply();

        Graphics.Blit(particlePositionTexture2D, particlePositionTextures[Read]);
        Destroy(particlePositionTexture2D);

        // Initialize velocities to zero
        ClearTexture(particleVelocityTextures[Read]);
    }

    private Vector3[] CreateParticlePositions()
    {
        Quaternion rotation = Quaternion.identity;
        Vector3 particleScale = new(particleRadius * 2, particleRadius * 2, particleRadius * 2);

        Vector3[] particlesPositions = new Vector3[desiredParticleCount];
        MeshProperties[] properties = new MeshProperties[desiredParticleCount];

        // Calculate dimensions for grid-based particle placement
        float particleSpacing = 2.0f * particleRadius;
        float particleVolume = particleMass * restDensity;

        // Calculate required volume based on particle count and spacing
        float totalVolume = desiredParticleCount * particleVolume;
        float sideLengthXYZ = Mathf.Pow(totalVolume, 1f / 3f);

        // Adjust simulation scale to fit particles
        transform.localScale = new Vector3(sideLengthXYZ, sideLengthXYZ / initialFluidHeight, sideLengthXYZ);

        // Calculate number of particles in each dimension
        int numX = Mathf.FloorToInt(transform.localScale.x / particleSpacing);
        int numY = Mathf.FloorToInt(transform.localScale.y * initialFluidHeight / particleSpacing);
        int numZ = Mathf.FloorToInt(transform.localScale.z / particleSpacing);

        // Offset to center particles in the volume
        Vector3 startPos = transform.position - transform.localScale / 2.0f;

        // Create particles in a grid pattern
        int particleIndex = 0;
        for (int y = 0; y < numY && particleIndex < desiredParticleCount; y++)
        {
            for (int x = 0; x < numX && particleIndex < desiredParticleCount; x++)
            {
                for (int z = 0; z < numZ && particleIndex < desiredParticleCount; z++)
                {
                    Vector3 position = startPos + new Vector3(
                        (x + 0.5f) * particleSpacing,
                        (y + 0.5f) * particleSpacing,
                        (z + 0.5f) * particleSpacing
                    );

                    MeshProperties props = new()
                    {
                        Mat = Matrix4x4.TRS(position, rotation, particleScale),
                        Color = Color.blue
                    };

                    particlesPositions[particleIndex] = position;
                    properties[particleIndex] = props;
                    particleIndex++;
                }
            }
        }

        _particleMeshPropertiesBuffer.SetData(properties);

        return particlesPositions;
    }

    private void InitializeBucketBuffer()
    {
        // Initialize bucket buffer with maximum possible particles per cell
        int totalBucketSize = bucketResolution * bucketResolution * bucketResolution * MaxParticlesPerVoxel;
        bucketBuffer = new ComputeBuffer(totalBucketSize, sizeof(uint));
    }

    public void InitializeCameraOrbit()
    {
        var cameraOrbit = Camera.main.GetComponent<CameraOrbit>();

        cameraOrbit.distance = transform.localScale.y + 10;
    }

    #endregion

    #region Update

    private void BucketGeneration()
    {
        // Clear bucket buffer with particle count as empty marker
        int totalBucketSize = bucketResolution * bucketResolution * bucketResolution * MaxParticlesPerVoxel;
        uint[] clearData = new uint[totalBucketSize];
        for (int i = 0; i < totalBucketSize; i++)
            clearData[i] = (uint)desiredParticleCount;
        bucketBuffer.SetData(clearData);

        // Set shader parameters
        bucketShader.SetBuffer(0, ShaderIDs.Bucket, bucketBuffer);
        bucketShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTextures[Read]);
        bucketShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        bucketShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleTextureResolution, particleTextureResolution));
        bucketShader.SetVector(ShaderIDs.BucketResolution, new(bucketWidth, bucketHeight, bucketDepth));
        bucketShader.SetVector(ShaderIDs.SimOrigin, transform.position - transform.localScale / 2);
        bucketShader.SetVector(ShaderIDs.SimScale, transform.localScale);

        // Calculate dispatch size for 2D thread groups
        int threadGroupsX = Mathf.CeilToInt((float)particleTextureResolution / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particleTextureResolution / NumThreads);
        bucketShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
    }

    private void DensityCalculation()
    {
        // Clear density texture
        ClearTexture(particleDensityTexture);

        // Set shader parameters
        densityShader.SetTexture(0, ShaderIDs.ParticleDensityTexture, particleDensityTexture);
        densityShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTextures[Read]);
        densityShader.SetBuffer(0, ShaderIDs.Bucket, bucketBuffer);

        densityShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        densityShader.SetFloat(ShaderIDs.ParticleMass, particleMass);
        densityShader.SetFloat(ShaderIDs.EffectiveRadius2, effectiveRadius * effectiveRadius);
        densityShader.SetFloat(ShaderIDs.EffectiveRadius9, Mathf.Pow(effectiveRadius, 9));
        densityShader.SetVector(ShaderIDs.SimOrigin, transform.position - transform.localScale / 2);
        densityShader.SetVector(ShaderIDs.SimScale, transform.localScale);
        densityShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleTextureResolution, particleTextureResolution));

        // Calculate dispatch size for 2D thread groups
        int threadGroupsX = Mathf.CeilToInt((float)particleTextureResolution / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particleTextureResolution / NumThreads);
        densityShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
    }

    private void UpdateVelocityAndPosition()
    {
        velPosShader.SetTexture(0, ShaderIDs.ParticlePositionTextureWrite, particlePositionTextures[Write]);
        velPosShader.SetTexture(0, ShaderIDs.ParticleVelocityTextureWrite, particleVelocityTextures[Write]);
        velPosShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTextures[Read]);
        velPosShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, particleVelocityTextures[Read]);
        velPosShader.SetTexture(0, ShaderIDs.ParticleDensityTexture, particleDensityTexture);
        velPosShader.SetBuffer(0, ShaderIDs.Bucket, bucketBuffer);
        velPosShader.SetBuffer(0, ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        velPosShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        velPosShader.SetFloat(ShaderIDs.EffectiveRadius, effectiveRadius);
        velPosShader.SetFloat(ShaderIDs.EffectiveRadius6, Mathf.Pow(effectiveRadius, 6));
        velPosShader.SetFloat(ShaderIDs.ParticleMass, particleMass);
        velPosShader.SetFloat(ShaderIDs.TimeStep, timeStep);
        velPosShader.SetFloat(ShaderIDs.Viscosity, viscosity);
        velPosShader.SetFloat(ShaderIDs.GasConst, gasConstant);
        velPosShader.SetFloat(ShaderIDs.RestDensity, restDensity);
        velPosShader.SetFloat(ShaderIDs.StiffnessCoeff, stiffnessCoeff);
        velPosShader.SetFloat(ShaderIDs.DampingCoeff, dampingCoeff);
        velPosShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleTextureResolution, particleTextureResolution));
        velPosShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        velPosShader.SetVector(ShaderIDs.SimOrigin, transform.position - transform.localScale / 2);
        velPosShader.SetVector(ShaderIDs.SimScale, transform.localScale);
        velPosShader.SetVector(ShaderIDs.ParticleScale, new(particleRadius, particleRadius, particleRadius));

        // Calculate dispatch size for 2D thread groups
        int threadGroupsX = Mathf.CeilToInt((float)particleTextureResolution / NumThreads);
        int threadGroupsY = Mathf.CeilToInt((float)particleTextureResolution / NumThreads);
        velPosShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);

        Swap(particlePositionTextures);
        Swap(particleVelocityTextures);
    }

    #endregion

    #region Texture Helpers

    private void Swap(RenderTexture[] textures)
    {
        (textures[Write], textures[Read]) = (textures[Read], textures[Write]);
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

    private void ClearTexture(RenderTexture texture)
    {
        if (texture.dimension == TextureDimension.Tex2D)
        {
            clearShader.SetTexture(0, ShaderIDs.Texture2D, texture);
            int threadGroupsX = Mathf.CeilToInt((float)texture.width / NumThreads);
            int threadGroupsY = Mathf.CeilToInt((float)texture.height / NumThreads);
            clearShader.Dispatch(0, threadGroupsX, threadGroupsY, 1);
        }
        else if (texture.dimension == TextureDimension.Tex3D)
        {
            clearShader.SetTexture(1, ShaderIDs.Texture3D, texture);
            int threadGroupsX = Mathf.CeilToInt((float)texture.width / NumThreads);
            int threadGroupsY = Mathf.CeilToInt((float)texture.height / NumThreads);
            int threadGroupsZ = Mathf.CeilToInt((float)texture.volumeDepth / NumThreads);
            clearShader.Dispatch(1, threadGroupsX, threadGroupsY, threadGroupsZ);
        }
    }

    #endregion
}
