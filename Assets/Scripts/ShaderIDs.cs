using UnityEngine;

public class ShaderIDs
{
    public static readonly int Texture2D = Shader.PropertyToID("_Texture2D");
    public static readonly int Texture3D = Shader.PropertyToID("_Texture3D");
    public static readonly int Properties = Shader.PropertyToID("_Properties");
    public static readonly int GridResolution = Shader.PropertyToID("_GridResolution");
    public static readonly int ParticleResolution = Shader.PropertyToID("_ParticleResolution");
    public static readonly int ParticlePositionTexture = Shader.PropertyToID("_ParticlePositionTexture");
    public static readonly int ParticleVelocityTexture = Shader.PropertyToID("_ParticleVelocityTexture");
    public static readonly int TimeStep = Shader.PropertyToID("_TimeStep");
    public static readonly int MeshTransform = Shader.PropertyToID("_MeshTransform");
    public static readonly int SimOrigin = Shader.PropertyToID("_SimOrigin");
    public static readonly int CellSize = Shader.PropertyToID("_CellSize");
    public static readonly int TriangleCount = Shader.PropertyToID("_TriangleCount");
    public static readonly int DistanceTexture = Shader.PropertyToID("_DistanceTexture");
    public static readonly int Vertices = Shader.PropertyToID("_Vertices");
    public static readonly int Triangles = Shader.PropertyToID("_Triangles");
    public static readonly int Bucket = Shader.PropertyToID("_Bucket");
    public static readonly int BucketResolution = Shader.PropertyToID("_BucketResolution");
    public static readonly int SimScale = Shader.PropertyToID("_SimScale");
    public static readonly int ParticleDensityTexture = Shader.PropertyToID("_ParticleDensityTexture");
    public static readonly int ParticleMass = Shader.PropertyToID("_ParticleMass");
    public static readonly int EffectiveRadius = Shader.PropertyToID("_EffectiveRadius");
    public static readonly int NumParticles = Shader.PropertyToID("_NumParticles");
}
