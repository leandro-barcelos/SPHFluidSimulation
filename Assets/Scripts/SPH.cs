using System;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;

public class SPH : MonoBehaviour
{
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

    [Header("Parameters")]
    [Range(0f, 1f / 60f)] public float timeStep = 1f / 60f;
    public int desiredParticleCount = 1000;
    [Range(0.01f, 1f)] public float spawnAreaFillRate;

    [Header("Rendering")]
    public float occlusionRange;
    public Material particleMaterial;
    public bool renderParticles;

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

    private void Start()
    {
        camera = Camera.main;
        cameraOrbit = camera.GetComponent<CameraOrbit>();
        InitializeBoxes();

        // Adjust particle resolution for better distribution
        var particlesWidth = Mathf.NextPowerOfTwo(Mathf.CeilToInt(Mathf.Sqrt(desiredParticleCount)));
        var particlesHeight = Mathf.CeilToInt((float)desiredParticleCount / particlesWidth);

        var particleCount = particlesWidth * particlesHeight;

        // Adjust particle radius based on grid resolution
        Vector3[] particlesPositions = CreateParticlePositions(particleCount);
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
        var particleRadius = Mathf.Pow(3 * particleVol / (4 * Mathf.PI), 1f / 3f) * spawnAreaFillRate;

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
    }

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
}
