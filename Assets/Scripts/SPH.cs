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

    [Header("Initialization")]
    [Range(1024, 4194304)] public int particleNumber = 1024;
    [Range(1, 256)] public int bucketResolution = 256;
    [Range(0.01f, 1f)] public float damFillRate = 0.5f;

    [Header("Parameters")]
    [Range(0.01f, 0.1f)] public float viscosity = 0.01f;
    [Range(1f, 2f)] public float restDensity = 1.5f;
    [Range(100f, 500f)] public float gasConstant = 150.0f;
    [Range(1000f, 10000f)] public float stiffnessCoeff = 5000.0f;
    [Range(1f, 50f)] public float dampingCoeff = 10.0f;

    [Header("Rendering")]
    public float occlusionRange;
    [Range(0.001f, 1f)] public float particleRadius;
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
    private Matrix4x4 simTRS;

    // Shaders
    private ComputeShader bucketShader, clearShader, densityShader, velPosShader, updateMeshPropertiesShader;
    private int threadGroups;

    #endregion

    #region Unity Functions

    private void Start()
    {
        particleNumber = Mathf.NextPowerOfTwo(particleNumber);
        particleTextureResolution = (int)Mathf.Sqrt(particleNumber);

        InitShaders();

        CreateParticleTextures();

        InitializeParticles();

        InitializeBucketBuffer();
    }

    private void Update()
    {
        simTRS = transform.localToWorldMatrix;

        BucketGeneration();
        DensityCalculation();

        for (var i = 0; i < 5; i++)
            UpdateVelocityAndPosition(Time.deltaTime / 25);

        UpdateMeshProperties();

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
        updateMeshPropertiesShader = Resources.Load<ComputeShader>("UpdateMeshProperties");

        threadGroups = Mathf.CeilToInt((float)particleTextureResolution / NumThreads);
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
        args[1] = (uint)particleNumber;
        args[2] = _particleMesh.GetIndexStart(0);
        args[3] = _particleMesh.GetBaseVertex(0);

        _particleArgsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _particleArgsBuffer.SetData(args);

        _particleMeshPropertiesBuffer = new ComputeBuffer(particleNumber, MeshProperties.Size());

        Vector3[] particlesPositions = CreateParticlePositions();

        particleMaterial.SetBuffer(ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        // Set texture to positions
        Texture2D particlePositionTexture2D = new(particleTextureResolution, particleTextureResolution, TextureFormat.RGBAFloat, false);
        Color[] particleColors = new Color[particleNumber];

        for (int i = 0; i < particleNumber; i++)
        {
            Vector3 pos = particlesPositions[i];
            particleColors[i] = new Color(pos.x, pos.y, pos.z, 1.0f);
        }

        particlePositionTexture2D.SetPixels(particleColors);
        particlePositionTexture2D.Apply();

        Graphics.Blit(particlePositionTexture2D, particlePositionTextures[Read]);
        Destroy(particlePositionTexture2D);

        // Initialize velocities to zero
        ClearTexture2D(particleVelocityTextures[Read]);
    }

    private Vector3[] CreateParticlePositions()
    {
        Quaternion rotation = Quaternion.identity;
        Vector3 particleScale = new(particleRadius * 2, particleRadius * 2, particleRadius * 2);

        Vector3[] particlesPositions = new Vector3[particleNumber];
        MeshProperties[] properties = new MeshProperties[particleNumber];

        var particlePerDim = Mathf.CeilToInt(Mathf.Pow(particleNumber / damFillRate, 1f / 3f));

        int xSize = Mathf.CeilToInt(particlePerDim * damFillRate);
        int ySize = particlePerDim;
        int zSize = particlePerDim;

        float particleCubeSize = 1f / particlePerDim;

        for (var i = 0; i < particleNumber; i++)
        {
            var position = new Vector3(
                particleCubeSize / 2f + i / (zSize * ySize) * damFillRate / xSize,
                particleCubeSize / 2f + i / zSize % ySize * 0.9f / ySize,
                particleCubeSize / 2f + i % zSize * 1f / zSize
            );

            position += new Vector3(
                Mathf.PerlinNoise(position.x, i) * 2 - 1,
                Mathf.PerlinNoise(position.y, i) * 2 - 1,
                Mathf.PerlinNoise(position.z, i) * 2 - 1
            ) * particleCubeSize;

            particlesPositions[i] = position;

            MeshProperties props = new()
            {
                Mat = Matrix4x4.TRS(position, rotation, particleScale),
                Color = Color.blue
            };

            properties[i] = props;
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

    #endregion

    #region Update

    private void BucketGeneration()
    {
        // Clear bucket buffer with particle count as empty marker
        int totalBucketSize = bucketResolution * bucketResolution * bucketResolution * MaxParticlesPerVoxel;
        uint[] clearData = new uint[totalBucketSize];
        for (int i = 0; i < totalBucketSize; i++)
            clearData[i] = (uint)particleNumber;
        bucketBuffer.SetData(clearData);

        // Set shader parameters
        bucketShader.SetBuffer(0, ShaderIDs.Bucket, bucketBuffer);
        bucketShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTextures[Read]);
        bucketShader.SetInt(ShaderIDs.NumParticles, particleNumber);
        bucketShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        bucketShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleTextureResolution, particleTextureResolution));

        bucketShader.Dispatch(0, threadGroups, threadGroups, 1);
    }

    private void DensityCalculation()
    {
        // Clear density texture
        ClearTexture2D(particleDensityTexture);

        // Set shader parameters
        densityShader.SetTexture(0, ShaderIDs.ParticleDensityTexture, particleDensityTexture);
        densityShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTextures[Read]);
        densityShader.SetBuffer(0, ShaderIDs.Bucket, bucketBuffer);

        densityShader.SetInt(ShaderIDs.NumParticles, particleNumber);
        densityShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        densityShader.SetFloat(ShaderIDs.ParticleMass, particleMass);
        densityShader.SetFloat(ShaderIDs.EffectiveRadius2, effectiveRadius * effectiveRadius);
        densityShader.SetFloat(ShaderIDs.EffectiveRadius9, Mathf.Pow(effectiveRadius, 9));
        densityShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleTextureResolution, particleTextureResolution));

        densityShader.Dispatch(0, threadGroups, threadGroups, 1);
    }

    private void UpdateVelocityAndPosition(float dt)
    {
        velPosShader.SetTexture(0, ShaderIDs.ParticlePositionTextureWrite, particlePositionTextures[Write]);
        velPosShader.SetTexture(0, ShaderIDs.ParticleVelocityTextureWrite, particleVelocityTextures[Write]);
        velPosShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTextures[Read]);
        velPosShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, particleVelocityTextures[Read]);
        velPosShader.SetTexture(0, ShaderIDs.ParticleDensityTexture, particleDensityTexture);
        velPosShader.SetBuffer(0, ShaderIDs.Bucket, bucketBuffer);

        velPosShader.SetInt(ShaderIDs.NumParticles, particleNumber);
        velPosShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        velPosShader.SetFloat(ShaderIDs.EffectiveRadius, effectiveRadius);
        velPosShader.SetFloat(ShaderIDs.EffectiveRadius6, Mathf.Pow(effectiveRadius, 6));
        velPosShader.SetFloat(ShaderIDs.ParticleMass, particleMass);
        velPosShader.SetFloat(ShaderIDs.TimeStep, dt);
        velPosShader.SetFloat(ShaderIDs.Viscosity, viscosity);
        velPosShader.SetFloat(ShaderIDs.GasConst, gasConstant);
        velPosShader.SetFloat(ShaderIDs.RestDensity, restDensity);
        velPosShader.SetFloat(ShaderIDs.StiffnessCoeff, stiffnessCoeff);
        velPosShader.SetFloat(ShaderIDs.DampingCoeff, dampingCoeff);
        velPosShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleTextureResolution, particleTextureResolution));

        velPosShader.Dispatch(0, threadGroups, threadGroups, 1);

        Swap(particlePositionTextures);
        Swap(particleVelocityTextures);
    }

    private void UpdateMeshProperties()
    {
        updateMeshPropertiesShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, particlePositionTextures[Read]);
        updateMeshPropertiesShader.SetTexture(0, ShaderIDs.ParticleDensityTexture, particleDensityTexture);
        updateMeshPropertiesShader.SetBuffer(0, ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        updateMeshPropertiesShader.SetFloat(ShaderIDs.RestDensity, restDensity);
        updateMeshPropertiesShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(particleTextureResolution, particleTextureResolution));
        updateMeshPropertiesShader.SetVector(ShaderIDs.ParticleScale, new(particleRadius, particleRadius, particleRadius));
        updateMeshPropertiesShader.SetMatrix(ShaderIDs.SimTRS, simTRS);

        updateMeshPropertiesShader.Dispatch(0, threadGroups, threadGroups, 1);
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

    private void ClearTexture2D(RenderTexture texture)
    {
        if (texture.dimension != TextureDimension.Tex2D)
            return;

            int threadGroupsX = Mathf.CeilToInt((float)texture.width / NumThreads);
            int threadGroupsY = Mathf.CeilToInt((float)texture.height / NumThreads);

        int kernel;

        switch (texture.format)
        {
            case RenderTextureFormat.RFloat:
                kernel = clearShader.FindKernel("ClearFloat");

                clearShader.SetTexture(kernel, ShaderIDs.Texture2D, texture);
                clearShader.Dispatch(kernel, threadGroupsX, threadGroupsY, 1);
                break;
            case RenderTextureFormat.ARGBFloat:
                kernel = clearShader.FindKernel("ClearFloat4");

                clearShader.SetTexture(kernel, ShaderIDs.Texture2D4, texture);
                clearShader.Dispatch(kernel, threadGroupsX, threadGroupsY, 1);
                break;
        }
    }

    #endregion
}
