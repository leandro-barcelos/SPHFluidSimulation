using UnityEngine;
using UnityEngine.Rendering;

public class SphFluidSimulation : MonoBehaviour
{
    #region Constants

    const int NumThreads = 32;
    const int MaxParticlesPerVoxel = 32;
    const int Read = 0;
    const int Write = 1;

    #endregion

    #region Auxiliary Structures
    private struct MeshProperties
    {
        // ReSharper disable once NotAccessedField.Local
        public Matrix4x4 Mat;
        // ReSharper disable once NotAccessedField.Local
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
    [Range(0, 2)] public int preset = 1;
    [Range(1024, 4194304)] public int particleNumber = 1024;
    [Range(1, 256)] public int bucketResolution = 256;
    [Range(0.01f, 1f)] public float damFillRate = 0.5f;

    [Header("Parameters")]
    [Range(0f, 0.1f)] public float viscosity = 0.01f;
    [Range(0f, 5f)] public float restDensity = 1.5f;
    [Range(1f, 5000f)] public float gasConstant = 150.0f;
    [Range(1000f, 10000f)] public float stiffnessCoefficient = 5000.0f;
    [Range(1f, 50f)] public float dampingCoefficient = 10.0f;

    [Header("Rendering")]
    public float occlusionRange;
    [Range(0.001f, 1f)] public float particleRadius;
    public Material particleMaterial;
    public bool renderParticles;
    [Range(0f, 1000f)] public float lowSpeed;
    [Range(0f, 1000f)] public float highSpeed;

    #endregion

    #region Private

    // Particle
    private RenderTexture[] _particlePositionTextures, _particleVelocityTextures;
    private RenderTexture _particleDensityTexture;
    private int _particleTextureResolution;
    private float _effectiveRadius;
    private float _particleMass;

    // Bucket
    private ComputeBuffer _bucketBuffer;

    // Rendering
    private Mesh _particleMesh;
    private ComputeBuffer _particleMeshPropertiesBuffer, _particleArgsBuffer;
    private Bounds _bounds;

    // Shaders
    private ComputeShader _bucketShader, _clearShader, _densityShader, _velPosShader, _updateMeshPropertiesShader, _initParticlesShader;
    private int _threadGroups;

    #endregion

    #region Unity Functions

    private void Start()
    {
        particleNumber = Mathf.NextPowerOfTwo(particleNumber);
        _particleTextureResolution = (int)Mathf.Sqrt(particleNumber);

        InitShaders();

        CreateParticleTextures();

        InitializeParticles();

        InitializeBucketBuffer();
    }

    private void Update()
    {
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
        _bucketBuffer?.Release();
        _particlePositionTextures[Read].Release();
        _particleVelocityTextures[Read].Release();
        _particlePositionTextures[Write].Release();
        _particleVelocityTextures[Write].Release();
        _particleDensityTexture?.Release();
    }

    #endregion

    #region Initializations

    private void InitShaders()
    {
        _clearShader = Resources.Load<ComputeShader>("Clear");
        _bucketShader = Resources.Load<ComputeShader>("Bucket");
        _densityShader = Resources.Load<ComputeShader>("Density");
        _velPosShader = Resources.Load<ComputeShader>("VelPos");
        _updateMeshPropertiesShader = Resources.Load<ComputeShader>("UpdateMeshProperties");
        _initParticlesShader = Resources.Load<ComputeShader>("InitParticles");

        _threadGroups = Mathf.CeilToInt((float)_particleTextureResolution / NumThreads);
    }

    private void CreateParticleTextures()
    {
        // Create particle position textures
        _particlePositionTextures = new RenderTexture[2];

        _particlePositionTextures[Read] = CreateRenderTexture2D(_particleTextureResolution, _particleTextureResolution, RenderTextureFormat.ARGBFloat);
        _particlePositionTextures[Write] = CreateRenderTexture2D(_particleTextureResolution, _particleTextureResolution, RenderTextureFormat.ARGBFloat);

        // Create particle velocity textures
        _particleVelocityTextures = new RenderTexture[2];

        _particleVelocityTextures[Read] = CreateRenderTexture2D(_particleTextureResolution, _particleTextureResolution, RenderTextureFormat.ARGBFloat);
        _particleVelocityTextures[Write] = CreateRenderTexture2D(_particleTextureResolution, _particleTextureResolution, RenderTextureFormat.ARGBFloat);

        // Create density texture
        _particleDensityTexture = CreateRenderTexture2D(_particleTextureResolution, _particleTextureResolution, RenderTextureFormat.RFloat);

    }

    private void InitializeParticles()
    {
        _effectiveRadius = 1f / (bucketResolution - 1);

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

        _particleMass = damFillRate / particleNumber;

        // Init  Positions
        _initParticlesShader.SetInt(ShaderIDs.ParticleResolution, _particleTextureResolution);
        _initParticlesShader.SetFloat(ShaderIDs.DamFillRate, damFillRate);

        _initParticlesShader.SetTexture(preset, ShaderIDs.ParticlePositionTexture, _particlePositionTextures[Read]);

        _initParticlesShader.Dispatch(preset, _threadGroups, _threadGroups, 1);

        particleMaterial.SetBuffer(ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        // Initialize velocities to zero
        ClearTexture2D(_particleVelocityTextures[Read]);
    }

    private void InitializeBucketBuffer()
    {
        // Initialize bucket buffer with maximum possible particles per cell
        var totalBucketSize = bucketResolution * bucketResolution * bucketResolution * MaxParticlesPerVoxel;
        _bucketBuffer = new ComputeBuffer(totalBucketSize, sizeof(uint));
    }

    #endregion

    #region Update

    private void BucketGeneration()
    {
        // Set shader parameters
        _bucketShader.SetInt(ShaderIDs.NumParticles, particleNumber);
        _bucketShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        _bucketShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(_particleTextureResolution, _particleTextureResolution));

        // Clear bucket buffer with particle count as empty marker
        _bucketShader.SetBuffer(1, ShaderIDs.Bucket, _bucketBuffer);
        var bucketThreadGroups = Mathf.CeilToInt((float)bucketResolution / 10);

        _bucketShader.Dispatch(1, bucketThreadGroups, bucketThreadGroups, bucketThreadGroups);

        // Generate Bucket
        _bucketShader.SetBuffer(0, ShaderIDs.Bucket, _bucketBuffer);
        _bucketShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, _particlePositionTextures[Read]);

        _bucketShader.Dispatch(0, _threadGroups, _threadGroups, 1);
    }

    private void DensityCalculation()
    {
        // Clear density texture
        ClearTexture2D(_particleDensityTexture);

        // Set shader parameters
        _densityShader.SetTexture(0, ShaderIDs.ParticleDensityTexture, _particleDensityTexture);
        _densityShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, _particlePositionTextures[Read]);
        _densityShader.SetBuffer(0, ShaderIDs.Bucket, _bucketBuffer);

        _densityShader.SetInt(ShaderIDs.NumParticles, particleNumber);
        _densityShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        _densityShader.SetFloat(ShaderIDs.ParticleMass, _particleMass);
        _densityShader.SetFloat(ShaderIDs.EffectiveRadius2, _effectiveRadius * _effectiveRadius);
        _densityShader.SetFloat(ShaderIDs.EffectiveRadius9, Mathf.Pow(_effectiveRadius, 9));
        _densityShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(_particleTextureResolution, _particleTextureResolution));

        _densityShader.Dispatch(0, _threadGroups, _threadGroups, 1);
    }

    private void UpdateVelocityAndPosition(float dt)
    {
        _velPosShader.SetTexture(0, ShaderIDs.ParticlePositionTextureWrite, _particlePositionTextures[Write]);
        _velPosShader.SetTexture(0, ShaderIDs.ParticleVelocityTextureWrite, _particleVelocityTextures[Write]);
        _velPosShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, _particlePositionTextures[Read]);
        _velPosShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, _particleVelocityTextures[Read]);
        _velPosShader.SetTexture(0, ShaderIDs.ParticleDensityTexture, _particleDensityTexture);
        _velPosShader.SetBuffer(0, ShaderIDs.Bucket, _bucketBuffer);

        _velPosShader.SetInt(ShaderIDs.NumParticles, particleNumber);
        _velPosShader.SetInt(ShaderIDs.BucketResolution, bucketResolution);
        _velPosShader.SetFloat(ShaderIDs.EffectiveRadius, _effectiveRadius);
        _velPosShader.SetFloat(ShaderIDs.EffectiveRadius6, Mathf.Pow(_effectiveRadius, 6));
        _velPosShader.SetFloat(ShaderIDs.ParticleMass, _particleMass);
        _velPosShader.SetFloat(ShaderIDs.TimeStep, dt);
        _velPosShader.SetFloat(ShaderIDs.Viscosity, viscosity);
        _velPosShader.SetFloat(ShaderIDs.GasConst, gasConstant);
        _velPosShader.SetFloat(ShaderIDs.RestDensity, restDensity);
        _velPosShader.SetFloat(ShaderIDs.StiffnessCoeff, stiffnessCoefficient);
        _velPosShader.SetFloat(ShaderIDs.DampingCoeff, dampingCoefficient);
        _velPosShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(_particleTextureResolution, _particleTextureResolution));

        _velPosShader.Dispatch(0, _threadGroups, _threadGroups, 1);

        Swap(_particlePositionTextures);
        Swap(_particleVelocityTextures);
    }

    private void UpdateMeshProperties()
    {
        _updateMeshPropertiesShader.SetTexture(0, ShaderIDs.ParticlePositionTexture, _particlePositionTextures[Read]);
        _updateMeshPropertiesShader.SetTexture(0, ShaderIDs.ParticleVelocityTexture, _particleVelocityTextures[Read]);
        _updateMeshPropertiesShader.SetBuffer(0, ShaderIDs.Properties, _particleMeshPropertiesBuffer);

        _updateMeshPropertiesShader.SetFloat(ShaderIDs.HighSpeed, highSpeed);
        _updateMeshPropertiesShader.SetFloat(ShaderIDs.LowSpeed, lowSpeed);
        _updateMeshPropertiesShader.SetVector(ShaderIDs.ParticleResolution, new Vector2(_particleTextureResolution, _particleTextureResolution));
        _updateMeshPropertiesShader.SetVector(ShaderIDs.ParticleScale, new Vector4(particleRadius, particleRadius, particleRadius));
        _updateMeshPropertiesShader.SetMatrix(ShaderIDs.SimTRS, transform.localToWorldMatrix);

        _updateMeshPropertiesShader.Dispatch(0, _threadGroups, _threadGroups, 1);
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
                kernel = _clearShader.FindKernel("ClearFloat");

                _clearShader.SetTexture(kernel, ShaderIDs.Texture2D, texture);
                _clearShader.Dispatch(kernel, threadGroupsX, threadGroupsY, 1);
                break;
            case RenderTextureFormat.ARGBFloat:
                kernel = _clearShader.FindKernel("ClearFloat4");

                _clearShader.SetTexture(kernel, ShaderIDs.Texture2D4, texture);
                _clearShader.Dispatch(kernel, threadGroupsX, threadGroupsY, 1);
                break;
        }
    }

    #endregion
}
