// https://www.alglib.net/translator/man/manual.csharp.html#sub_spdmatrixcholesky

using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// AoS
public struct ClothNode
{
    public float x, y, z;
    public float vx, vy, vz;
    public int ax, ay, az;
    public float mass;
}

// SoA
public struct ClothNodesAttrArray
{
    public Vector3 [,] pos;
    public Vector3 [,] vel;
    public Vector3 [,] acc;
}

public class OptMethod : MonoBehaviour
{
    public GameObject hairPrefab, colliderPrefab;
    public LampOptMethod theLamp;
    public Texture2D weightMap;

    static public int simulationSteps;
    public static int nHairs, nNodesPerHair;

    // array of ClothNode Struts, Object Oriented design.
    static ClothNode[] clothNodesArray;
    static GameObject[] clothVertArray;

    // ClothNodesAttrArray
    ClothNodesAttrArray clothAttrsArray; 

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

    private void Awake()
    {
        initData();
    }

    // Start is called before the first frame update
    void Start()
    {
        initGeo();
    }

    // Update is called once per frame
    void Update()
    {
        tickSimulation(1);
        updateClothVertexPositions();
    }

    void initData()
    {
        simulationSteps = 40;

        // hair nodes
        nHairs = 32;
        nNodesPerHair = 32;

        // cloth Node
        clothNodesArray = new ClothNode[nHairs * nNodesPerHair];

        //cloth NodeList
        clothAttrsArray.pos = new Vector3 [nHairs, nNodesPerHair];
        clothAttrsArray.vel = new Vector3 [nHairs, nNodesPerHair];
        clothAttrsArray.acc = new Vector3 [nHairs, nNodesPerHair];

        // simulation variables
        nodeDistance = 0.5f;    // Initial Node distance apart.
        dPosition = 0.0004f;    // Euler method integration ratio for speed.
        dVelocity = 1.0f;       // Euler method integration ratio for acceleration.
        forceDecay = 0.0000f;
        velocityDecay = 0.999f;
        gravity = 0.1f;
        stiffness = 6.0f;
        maxTravelDistance = 5.0f;
        bendingStiffness = 0.1f;

        // an Array of Node struct
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                // AoS
                int nodeIndex = i * nNodesPerHair + j;
                clothNodesArray[nodeIndex].x = nodeDistance * (i - nHairs / 2);
                clothNodesArray[nodeIndex].y = -nodeDistance * (j - nNodesPerHair / 2);
                clothNodesArray[nodeIndex].z = 0.0f;
                clothNodesArray[nodeIndex].vx = 0.0f;
                clothNodesArray[nodeIndex].vy = 0.0f;
                clothNodesArray[nodeIndex].vz = 0.0f;
                clothNodesArray[nodeIndex].ax = 0;
                clothNodesArray[nodeIndex].ay = 0;
                clothNodesArray[nodeIndex].az = 0;

                // SoA
                clothAttrsArray.pos[i, j].x = nodeDistance * (i - nHairs / 2);
                clothAttrsArray.pos[i, j].y = -nodeDistance * (j - nNodesPerHair / 2);
                clothAttrsArray.pos[i, j].z = 0.0f;
                clothAttrsArray.vel[i, j].x = 0.0f;
                clothAttrsArray.vel[i, j].y = 0.0f;
                clothAttrsArray.vel[i, j].z = 0.0f;
                clothAttrsArray.acc[i, j].x = 0.0f;
                clothAttrsArray.acc[i, j].y = 0.0f;
                clothAttrsArray.acc[i, j].z = 0.0f;

                // sample the weight from the weight map
                float u = i * 1.0f / (nHairs);
                float v = 1.0f - j * 1.0f / (nNodesPerHair);
                clothNodesArray[nodeIndex].mass = 1.0f + weightMap.GetPixelBilinear(u, v).grayscale;
            }
        }
    }

    void initGeo()
    {
        // instantiate hair objects
        clothVertArray = new GameObject[nHairs * nNodesPerHair];
        for (int i = 0; i < clothVertArray.Length; i++)
        {
            Vector3 location = new Vector3(clothNodesArray[i].x, clothNodesArray[i].y, clothNodesArray[i].z);
            var newitem = Instantiate(hairPrefab, location, Quaternion.identity);
            clothVertArray[i] = newitem;
        }
    }

    void updateClothVertexPositions()
    {
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++) 
            {
                clothVertArray[i*nNodesPerHair + j].transform.position = clothAttrsArray.pos[i, j];
            }
        }
    }

    void tickSimulation(int simulationIndex)
    {
        switch (simulationIndex)
        {
            case 1:
                explicitEulerSimulationMethod();
                // p [matrix of vector3 # of nodes I have] = p+v*h
                // v = v+ a*m
                // for every node. 
                break;
            case 2:
                localGlobalSimulationMethod();
                // finish implementation? understand what is going on here
                // Expectation: continue this research. 
                // by early fall, be able 
                
                break;
            case 3:
                GPUSimulationMethod();
                break;
            default:
                Debug.Log("Unknown Simulation Method!");
                break;
        }
    }

    void localGlobalSimulationMethod()
    {
        return;
    }

    void explicitEulerSimulationMethod()
    {
        for (int i = 0; i < simulationSteps; i++)
        {
            explicitEuler_tickVel();
            explicitEuler_tickForce();
            explicitEuler_tickIntegration();
        }
        return;
    }

    void GPUSimulationMethod()
    {
        return;
    }

    // kernel Force in CPU implementation 
    void explicitEuler_tickForce()
    {
        // calculate F among the nodes
        for (int i = 0; i <  nHairs; i++) {
            for (int j = 0; j < nNodesPerHair - 1; j++)
            {
                // spring force
                Vector3 currNext = new Vector3( clothAttrsArray.pos[i, j + 1].x - clothAttrsArray.pos[i, j].x,
                                                clothAttrsArray.pos[i, j + 1].y - clothAttrsArray.pos[i, j].y,
                                                clothAttrsArray.pos[i, j + 1].z - clothAttrsArray.pos[i, j].z);

                float dX = Mathf.Clamp(nodeDistance - currNext.magnitude, -maxTravelDistance, maxTravelDistance);
                Vector3 springForce = -stiffness * currNext.normalized * dX;

                // bending force
                Vector3 bendingForce = new Vector3(0.0f, 0.0f, 0.0f);

                if (j != 0) {
                    Vector3 currPrev = new Vector3( clothAttrsArray.pos[i, j - 1].x - clothAttrsArray.pos[i, j].x,
                                                    clothAttrsArray.pos[i, j - 1].y - clothAttrsArray.pos[i, j].y,
                                                    clothAttrsArray.pos[i, j - 1].z - clothAttrsArray.pos[i, j].z);

                    bendingForce = bendingStiffness * (currNext + currPrev);
                    clothAttrsArray.acc[i, j - 1] -= 0.5f * bendingForce;
                }

                clothAttrsArray.acc[i, j] += bendingForce;
                clothAttrsArray.acc[i, j] += springForce;
                clothAttrsArray.acc[i, j].y -= gravity;

                clothAttrsArray.acc[i, j + 1] -= springForce;
                clothAttrsArray.acc[i, j + 1] -= 0.5f * bendingForce ;
            }
        }
    }

    void explicitEuler_tickVel()
    {
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair - 1; j++)
            {
                Vector3 relativePos = new Vector3(  clothAttrsArray.pos[i, j + 1].x - clothAttrsArray.pos[i, j].x,
                                                    clothAttrsArray.pos[i, j + 1].y - clothAttrsArray.pos[i, j].y,
                                                    clothAttrsArray.pos[i, j + 1].z - clothAttrsArray.pos[i, j].z);

                Vector3 relativeVel = new Vector3(  clothAttrsArray.vel[i, j + 1].x - clothAttrsArray.vel[i, j].x,
                                                    clothAttrsArray.vel[i, j + 1].y - clothAttrsArray.vel[i, j].y,
                                                    clothAttrsArray.vel[i, j + 1].z - clothAttrsArray.vel[i, j].z);

                Vector3 tangentVel = Vector3.Dot(relativeVel, relativePos.normalized) * relativePos.normalized;
                Vector3 verticalVel = relativeVel - tangentVel;
                Vector3 exVel = 0.005f * tangentVel + 0.001f * verticalVel;

                clothAttrsArray.acc[i, j] += exVel;
                clothAttrsArray.acc[i, j + 1] -= exVel;
            }
        }
    }

    void explicitEuler_tickIntegration()
    {
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                if (j == 0)
                {
                    clothAttrsArray.pos[i, j].x = nodeDistance * (i - nHairs / 2);
                    clothAttrsArray.pos[i, j].y = -nodeDistance * (j - nNodesPerHair / 2);
                    clothAttrsArray.pos[i, j].z = 0.0f;
                    clothAttrsArray.vel[i, j].x = 0.0f;
                    clothAttrsArray.vel[i, j].y = 0.0f;
                    clothAttrsArray.vel[i, j].z = 0.0f;
                    clothAttrsArray.acc[i, j].x = 0.0f;
                    clothAttrsArray.acc[i, j].y = 0.0f;
                    clothAttrsArray.acc[i, j].z = 0.0f;
                }

                // Euler integration
                // v = v + ah
                clothAttrsArray.vel[i, j] += dVelocity * clothAttrsArray.acc[i, j];

                // p = v + vh
                clothAttrsArray.pos[i, j] += dPosition * clothAttrsArray.vel[i, j];

                clothAttrsArray.acc[i, j] *= forceDecay;
                clothAttrsArray.vel[i, j] *= velocityDecay;
            }
        }
    }
    void explicitEuler_tickLight()
    {

    }

    public void updateLightSettings(Vector3 lightHeadPos, Vector3 lightHeadLootAt)
    {
        headPos = lightHeadPos;
        headLookAtPos = lightHeadLootAt;
    }
}
