using UnityEngine;

public class ShaderIDs
{
    public static readonly int Texture2D = Shader.PropertyToID("_Texture2D");
    public static readonly int Texture2D4 = Shader.PropertyToID("_Texture2D4");
    public static readonly int Properties = Shader.PropertyToID("_Properties");
    public static readonly int ParticleResolution = Shader.PropertyToID("_ParticleResolution");
    public static readonly int ParticlePositionTexture = Shader.PropertyToID("_ParticlePositionTexture");
    public static readonly int ParticleVelocityTexture = Shader.PropertyToID("_ParticleVelocityTexture");
    public static readonly int TimeStep = Shader.PropertyToID("_TimeStep");
    public static readonly int Bucket = Shader.PropertyToID("_Bucket");
    public static readonly int BucketResolution = Shader.PropertyToID("_BucketResolution");
    public static readonly int ParticleDensityTexture = Shader.PropertyToID("_ParticleDensityTexture");
    public static readonly int ParticleMass = Shader.PropertyToID("_ParticleMass");
    public static readonly int EffectiveRadius = Shader.PropertyToID("_EffectiveRadius");
    public static readonly int EffectiveRadius2 = Shader.PropertyToID("_EffectiveRadius2");
    public static readonly int EffectiveRadius9 = Shader.PropertyToID("_EffectiveRadius9");
    public static readonly int EffectiveRadius6 = Shader.PropertyToID("_EffectiveRadius6");
    public static readonly int Viscosity = Shader.PropertyToID("_Viscosity");
    public static readonly int RestDensity = Shader.PropertyToID("_RestDensity");
    public static readonly int GasConst = Shader.PropertyToID("_GasConst");
    public static readonly int ParticleScale = Shader.PropertyToID("_ParticleScale");
    public static readonly int ParticlePositionTextureWrite = Shader.PropertyToID("_ParticlePositionTextureWrite");
    public static readonly int ParticleVelocityTextureWrite = Shader.PropertyToID("_ParticleVelocityTextureWrite");
    public static readonly int NumParticles = Shader.PropertyToID("_NumParticles");
    public static readonly int StiffnessCoeff = Shader.PropertyToID("_StiffnessCoeff");
    public static readonly int DampingCoeff = Shader.PropertyToID("_DampingCoeff");
    public static readonly int SimTRS = Shader.PropertyToID("_SimTRS");
    public static readonly int HighSpeed = Shader.PropertyToID("_HighSpeed");
    public static readonly int LowSpeed = Shader.PropertyToID("_LowSpeed");
    public static readonly int DamFillRate = Shader.PropertyToID("_DamFillRate");
}
