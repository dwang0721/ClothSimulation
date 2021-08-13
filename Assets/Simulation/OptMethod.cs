// Math Library
// https://numerics.mathdotnet.com/

// Math Plugin: https://github.com/GlitchEnzo/NuGetForUnity, following the steps to install MathNet:
// 1. add the NugetForUnity.unitypackage to your Unity Project.
// 2. Open the NuGet > Manage NuGet Packages Window.
// 3. Search the MathNet.Numerics and install.

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

// GPU method Data struct
public struct ClothNode
{
    public float x, y, z;
    public float vx, vy, vz;
    public int ax, ay, az;
    public float mass;
}

// For CPU method
public struct ClothNodesAttrArray
{
    public Vector3 [,] pos;
    public Vector3 [,] vel;
    public Vector3 [,] acc;
}

// For Matrix method
public struct ClothNodeMatrix
{
    Matrix<double> x, y, z;
    Matrix<double> vx, vy, vz;
    Matrix<double> ax, ay, az;
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
        //explicitMatrixSimulationMethod();
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
                // AoS GPU method Data struct
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

                // SoA CPU method Data struct
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
            // do it in dynaic way, so that we don't have to thru all iterations to reach destinations.
            explicitEuler_tickVel();
            explicitEuler_tickForce();
            explicitEuler_tickLight();
            explicitEuler_tickIntegration();
            // difference is small quit iteration. adjustment is small, quit. 
            // skip nodes. 
            // check if it should continue. 
            // define what is close?
            // do it linear 
            // integration may oscillate. Detect error/ oscillate, in the code, 
            // analysis, understand the relation to localGlobalSimulationMethod(); 
        }
    }

    void explicitMatrixSimulationMethod()
    {
        Matrix<double> A = DenseMatrix.OfArray(new double[,] {  {1,1,1,1},
                                                                {1,2,3,4},
                                                                {4,3,2,1}   });

        Matrix<double> B = DenseMatrix.OfArray(new double[,] {  {-1,-1,-1,-1},
                                                                {1,2,3,4},
                                                                {-4,-3,-2,-1}   });

        Debug.Log(A + B);

        //Vector3[,] c = new Vector3[,] {  {new Vector3(1,1,1), new Vector3(2,2,2)},
        //                                 {new Vector3(1,1,1), new Vector3(2,2,2)} };
        //Matrix<Vector3> C = c;
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
                Vector3 currNext = clothAttrsArray.pos[i, j + 1] - clothAttrsArray.pos[i, j];
                float dX = Mathf.Clamp(nodeDistance - currNext.magnitude, -maxTravelDistance, maxTravelDistance);
                Vector3 springForce = -stiffness * currNext.normalized * dX;

                // bending force
                Vector3 bendingForce = new Vector3(0.0f, 0.0f, 0.0f);

                if (j != 0) {
                    Vector3 currPrev = clothAttrsArray.pos[i, j - 1] - clothAttrsArray.pos[i, j];
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
                Vector3 relativePos = clothAttrsArray.pos[i, j + 1] - clothAttrsArray.pos[i, j];
                Vector3 relativeVel = clothAttrsArray.vel[i, j + 1] - clothAttrsArray.vel[i, j];
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
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                Vector3 lookAtDir = (headLookAtPos - headPos).normalized;
                Vector3 lightDir = (clothAttrsArray.pos[i, j] - headPos).normalized;

                float dotValue = Vector3.Dot(lookAtDir, lightDir);
                Vector3 lightForceVector = lightDir * lightForce;

                if (dotValue > Mathf.Cos(lightAngle * 0.0174533f * 0.5f))
                {
                    clothAttrsArray.acc[i, j] += lightForceVector; 
                }
            }
        }
    }

    public void updateLightLooAtSettings(Vector3 lightHeadPos, Vector3 lightHeadLootAt)
    {
        headPos = lightHeadPos;
        headLookAtPos = lightHeadLootAt;
    }

    public void updatelightForce(float f)
    {
        //shader.SetFloat("lightForce", lightForce);
        lightForce = f;
    }

    public void updatelightAngle(float a)
    {
        //shader.SetFloat("lightAngle", lightAngle);
        lightAngle = a;
    }
}
