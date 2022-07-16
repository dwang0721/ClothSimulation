using System.Collections;
using System.Collections.Generic;
using UnityEngine.EventSystems;
using UnityEngine;

public struct HairNode3D
{
    public float x, y, z;
    public float vx, vy, vz;
    public int ax, ay, az;
    public float mass; 
    public float nx, ny, nz;
    public int isPinned, dummy2, dummy3;
}

public struct ColliderNode3D
{
    public float x, y, z, r;
    public int ax, ay, az;
    public int dummy1;
}

public class CPU3D: MonoBehaviour
{
    public ComputeShader shader;    
    public GameObject hairPrefab, sphereColliderPrefab, planeColliderPrefab;
    public Lamp theLamp;
    public Texture2D weightMap;
    public int resolution;
    public bool isFreefallMode;
    public ClothMeshGenerator clothMesh;
    public float clothDebugNodeSize;

    static public int simulationSteps;

    // hair nodes,  hair is a vertical line of the cloth
    public static int nHairs, nNodesPerHair;
    static HairNode3D[] hairNodesArray;    
    static GameObject[] hairGeos;

    // sphere colliders
    static public int nColliders;
    static public float colliderRadius;
    static ColliderNode3D[] colliderNodeArrays;
    static public GameObject[] colliderGeos;

    // plane collider
    static public float planeColliderPosY = 0;

    // simulation parameters
    public static float nodeDistance;
    static float dPosition;             // speed ratio for euler's method
    static float dVelocity;             // force ratio for euler's method
    public static float forceDecay;
    public static float velocityDecay;
    public static float gravity;
    public static float stiffness;             // used for Hooke's law spring coefficient.
    public static float maxTravelDistance;     // the maximum distance node apart.
    public static float bendingStiffness;      // coefficient for bending forces

    public static float lightForce;
    public static float lightAngle;
    public static Vector3 headPos;
    public static Vector3 headLookAtPos;

    public static int resetPinnedNodeDistanceFlag;

    // kernel
    int forceKernel, lightKernel, velocityKernel, collisionKernel, eulerKernel;

    // buffer is an array
    ComputeBuffer hairNodeBuffer, colliderBuffer;

    void Start()
    {
        Debug.Assert(shader);
        Debug.Assert(hairPrefab);
        Debug.Assert(sphereColliderPrefab);
        Debug.Assert(planeColliderPrefab);
        Debug.Assert(theLamp);
        Debug.Assert(weightMap);
        Debug.Assert(resolution > 0);
        Debug.Assert(clothDebugNodeSize > 0);
        Debug.Assert(clothMesh);

        initData();
        initGeo();
        initBuffer();
        initShader();
        initClothMesh();
    }

    // Update is called once per frame
    void Update()
    {
        simulationOnGPU();
        updateHairGeoPositions();
        updateClothMesh();
        updateDataFromCollider();
    }

    void initData()
    {
        // kernel data
        simulationSteps = 40;

        // hair data
        nHairs = resolution;
        nNodesPerHair = resolution;
        hairNodesArray = new HairNode3D[ nHairs * nNodesPerHair ];

        // simulation variables
        nodeDistance = 0.49f;    // Initial Node distance apart.
        dPosition = 0.0004f;    // Euler method integration ratio for speed.
        dVelocity = 1.0f;       // Euler method integration ratio for acceleration.
        forceDecay = 0.0000f;
        velocityDecay = 0.999f; 
        gravity = 0.1f;
        stiffness = 6.0f;
        maxTravelDistance = 5.0f;
        bendingStiffness = 0.1f;

        // light force variables
        lightForce = 0.1f;
        lightAngle = 10.0f;
        headPos = theLamp.head.transform.position;
        headLookAtPos = theLamp.headLookAt;

        // a flag deciding if the node distance to affect the pinned node
        resetPinnedNodeDistanceFlag = 0;

        // nodes data
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                int nodeIndex = i * nNodesPerHair + j;

                hairNodesArray[nodeIndex].isPinned = 0;

                if (isFreefallMode)
                {
                    hairNodesArray[nodeIndex].x = nodeDistance * (i - nHairs / 2);
                    hairNodesArray[nodeIndex].y = 0.0f;
                    hairNodesArray[nodeIndex].z = -nodeDistance * (j - nNodesPerHair / 2);
                }
                else {
                    hairNodesArray[nodeIndex].x = nodeDistance * (i - nHairs / 2);
                    hairNodesArray[nodeIndex].y = -nodeDistance * (j - nNodesPerHair / 2);
                    hairNodesArray[nodeIndex].z = 0.0f;

                    // first row is intialized as pinned node
                    if (j == 0)
                    {
                        hairNodesArray[nodeIndex].isPinned = 1;
                    }
                }

                hairNodesArray[nodeIndex].vx = 0.0f;
                hairNodesArray[nodeIndex].vy = 0.0f;
                hairNodesArray[nodeIndex].vz = 0.0f;
                hairNodesArray[nodeIndex].ax = 0;
                hairNodesArray[nodeIndex].ay = 0;
                hairNodesArray[nodeIndex].az = 0;

                // sample the weight from the weight map
                float u = i * 1.0f / (nHairs);
                float v = 1.0f - j * 1.0f / (nNodesPerHair);                
                hairNodesArray[nodeIndex].mass = 1.0f + weightMap.GetPixelBilinear(u, v).grayscale;

                // dummy node, do nothing here
                hairNodesArray[nodeIndex].dummy2 = 0;
                hairNodesArray[nodeIndex].dummy3 = 0;
            }                
        }

        // plane collider data
        planeColliderPosY = -1 * resolution * nodeDistance / 2.0f;

        // sphere collider data
        nColliders = 3;
        colliderRadius = 2.5f;
        colliderNodeArrays = new ColliderNode3D[nColliders];
        for (int i = 0; i < nColliders; i++)
        {
            if (isFreefallMode)
            {
                colliderNodeArrays[i].x = colliderRadius * Mathf.Cos(i * Mathf.PI * 2.0f / nColliders );
                colliderNodeArrays[i].y = planeColliderPosY;
                colliderNodeArrays[i].z = colliderRadius * Mathf.Sin(i * Mathf.PI * 2.0f / nColliders );
                colliderNodeArrays[i].r = colliderRadius * Random.Range(1.0f, 2.0f);
            }
            else {
                colliderNodeArrays[i].x = colliderRadius * (i - nColliders / 2);
                colliderNodeArrays[i].y = -20.0f;
                colliderNodeArrays[i].z = 0.0f;
                colliderNodeArrays[i].r = colliderRadius;
            }
            
            colliderNodeArrays[i].ax = 0;
            colliderNodeArrays[i].ay = 0;
            colliderNodeArrays[i].az = 0;
        }
    }

    void initGeo()
    {
        // instantiate hair objects
        hairGeos = new GameObject[nHairs * nNodesPerHair];
        for (int i = 0; i < hairGeos.Length; i++)
        {
            Vector3 location = new Vector3(hairNodesArray[i].x, hairNodesArray[i].y, hairNodesArray[i].z);
            var newitem = Instantiate(hairPrefab, location, Quaternion.identity);
            newitem.transform.localScale = new Vector3(1, 1, 1) * clothDebugNodeSize;       
            newitem.GetComponent<NodeController>().setNodeIndex(i);
            newitem.GetComponent<NodeController>().setIsPinned(hairNodesArray[i].isPinned);
            hairGeos[i] = newitem;
        }

        //  instantiate sphere colliders
        colliderGeos = new GameObject[nColliders];
        for (int i = 0; i < nColliders; i++)
        {
            Vector3 location = new Vector3(colliderNodeArrays[i].x, colliderNodeArrays[i].y, colliderNodeArrays[i].z); 
            var newitem = Instantiate(sphereColliderPrefab, location, Quaternion.identity);
            newitem.transform.localScale = new Vector3(1, 1, 1) * 2 * colliderNodeArrays[i].r;
            colliderGeos[i] = newitem;
        }

        // instantiate plane collider
        Vector3 planeColliderLocation = new Vector3(0, planeColliderPosY, 0);
        var planeCollider = Instantiate(planeColliderPrefab, planeColliderLocation, Quaternion.identity);
        planeCollider.transform.localScale = new Vector3(1, 1, 1) * 4;
        planeCollider.GetComponent<Renderer>().material.color = new Color(1.0f, 1.0f, 1.0f, 0.5f);
    }

    void initClothMesh()
    {
        clothMesh.initMeshExplict(hairNodesArray, nHairs, nNodesPerHair);
    }

    void initBuffer()
    {
        hairNodeBuffer = new ComputeBuffer(hairNodesArray.Length, 4*16); // 4 byte for float or int, 10 items
        hairNodeBuffer.SetData(hairNodesArray);
        colliderBuffer = new ComputeBuffer(nColliders, 4*8);
        colliderBuffer.SetData(colliderNodeArrays);
    }

    public GameObject getCollider() {
        return colliderGeos.Length > 0 ? colliderGeos[0] : null;
    }

    public void updateNodeDistance(float dis) {
        shader.SetFloat("nodeDistance", dis);
    }

    public void updateStiffness(float stiff)
    {
        shader.SetFloat("stiffness", stiff);
    }

    public void updateGravity(float grav)
    {
        shader.SetFloat("gravity", grav);
    }

    public void updateMaxTravel(float maxTrav)
    {
        shader.SetFloat("maxTravelDistance", maxTrav);
    }

    public void updateBendStiff(float bendStiff)
    {
        shader.SetFloat("bendingStiffness", bendStiff);
    }

    public void updateFriction(float friction)
    {
        shader.SetFloat("velocityDecay", friction);
    }

    public void updateLightDirection(Vector3 headPos, Vector3 headLookAtPos)
    {
        shader.SetVector("headPos", headPos);
        shader.SetVector("headLookAtPos", headLookAtPos);
    }

    public void updatelightForce(float lightForce)
    {
        shader.SetFloat("lightForce", lightForce);
    }

    public void updatelightAngle(float lightAngle)
    {
        shader.SetFloat("lightAngle", lightAngle);
    }

    public void setHairNodesPinStatusAtIndex(int index, int pinStatus)
    {
        hairNodesArray[index].isPinned = pinStatus;
    }

    public void resetPinnedNodeDistance(int flag)
    {
        shader.SetInt("resetPinnedNodeDistanceFlag", flag);
        resetPinnedNodeDistanceFlag = flag;
    }

    void initShader() 
    {
        shader.SetInt("nNodesPerHair", nNodesPerHair);
        shader.SetInt("nHairs", nHairs);
        shader.SetInt("nColliders", nColliders);
        shader.SetFloat("nodeDistance", nodeDistance);
        shader.SetFloat("dPosition", dPosition);
        shader.SetFloat("dVelocity", dVelocity);
        shader.SetFloat("forceDecay", forceDecay);
        shader.SetFloat("velocityDecay", velocityDecay);
        shader.SetFloat("gravity", gravity);
        shader.SetFloat("stiffness", stiffness); 
        shader.SetFloat("maxTravelDistance", maxTravelDistance);
        shader.SetFloat("bendingStiffness", bendingStiffness);
        shader.SetFloat("planeColliderPosY", planeColliderPosY);

        shader.SetFloat("lightForce", lightForce);
        shader.SetFloat("lightAngle", lightAngle);
        shader.SetVector("headPos", headPos);
        shader.SetVector("headLookAtPos", headLookAtPos);

        shader.SetInt("floatToInt", 2 << 17);
        shader.SetFloat("intToFloat", 1f / (2 << 17));
        shader.SetInt("resetPinnedNodeDistance", 0);

        velocityKernel = shader.FindKernel("VelocityKernel");
        shader.SetBuffer(velocityKernel, "hairNodeBuffer", hairNodeBuffer);

        forceKernel = shader.FindKernel("ForceKernel");
        shader.SetBuffer(forceKernel, "hairNodeBuffer", hairNodeBuffer);

        lightKernel = shader.FindKernel("LightKernel");
        shader.SetBuffer(lightKernel, "hairNodeBuffer", hairNodeBuffer);

        collisionKernel = shader.FindKernel("CollisionKernel");
        shader.SetBuffer(collisionKernel, "hairNodeBuffer", hairNodeBuffer);
        shader.SetBuffer(collisionKernel, "colliderBuffer", colliderBuffer);

        eulerKernel = shader.FindKernel("EulerKernel");
        shader.SetBuffer(eulerKernel, "hairNodeBuffer", hairNodeBuffer);
    }

    void simulationOnGPU() 
    {
        int nThreadGrpsX = 2;
        int nThreadGrpsY = 2;

        // set node data
        hairNodeBuffer.SetData(hairNodesArray);

        // set buffer data from collider
        colliderBuffer.SetData(colliderNodeArrays);

        // update node positions in GPU
        for (int i = 0; i < simulationSteps; i++) // simulationSteps
        {
            shader.Dispatch(velocityKernel, nThreadGrpsX, nThreadGrpsY, 1); // exchange velocity among nodes
            shader.Dispatch(forceKernel, nThreadGrpsX, nThreadGrpsY, 1);    // update force from position ( hooke's law )
            shader.Dispatch(lightKernel, nThreadGrpsX, nThreadGrpsY, 1);    // light force
            shader.Dispatch(collisionKernel, nThreadGrpsX, nThreadGrpsY, 1);// update position, velocity and force upon collision
            shader.Dispatch(eulerKernel, nThreadGrpsX, nThreadGrpsY, 1);    // Euler's method to accumulate to new position.
        }

        if (resetPinnedNodeDistanceFlag == 1)
        {
            // always reset the flag to zero,
            // because we don't want to update pinned node unless the user clicks the reset button.
            resetPinnedNodeDistance(0); 
        }

        // set data to collider
        colliderBuffer.GetData(colliderNodeArrays);

        // set data to hairNode
        hairNodeBuffer.GetData(hairNodesArray);
    }

    void updateHairGeoPositions()
    {
        for (int i = 0; i < hairGeos.Length; i++)
        {
            Vector3 location = new Vector3(hairNodesArray[i].x, hairNodesArray[i].y, hairNodesArray[i].z);
            hairGeos[i].transform.position = location;
        }
    }

    void updateClothMesh()
    {
        clothMesh.updateMeshExplict(hairNodesArray);
    }

    void updateDataFromCollider()
    {
        for (int i = 0; i < nColliders; i++)
        {
            colliderNodeArrays[i].x = colliderGeos[i].transform.position.x;
            colliderNodeArrays[i].y = colliderGeos[i].transform.position.y;
            colliderNodeArrays[i].z = colliderGeos[i].transform.position.z;

            // don't move collider for now
            colliderNodeArrays[i].ax = 0;
            colliderNodeArrays[i].ay = 0;
            colliderNodeArrays[i].az = 0;

            // assuming the sphere scale is double of the radius.
            colliderNodeArrays[i].r = colliderGeos[i].transform.localScale.x / 2.0f;
        }
    }

    void OnDestroy()
    {                  
        hairNodeBuffer.Release();
        colliderBuffer.Release();
    }
}
