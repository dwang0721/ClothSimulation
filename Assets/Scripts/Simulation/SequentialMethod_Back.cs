
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

// For CPU method
public struct ClothNodesAttrArray
{
    public Vector3[,] pos;
    public Vector3[,] vel;
    public Vector3[,] acc;
}

public class SequentialMethod : MonoBehaviour
{
    // stuff can be seen from the UI
    public GameObject hairPrefab, colliderPrefab;
    public LampOptMethod theLamp;
    public Texture2D weightMap;

    static GameObject[] clothVertArray;

    // CPU Method Data Structure.
    ClothNodesAttrArray clothAttrsArray;

    static public int simulationSteps;
    static public int nHairs, nNodesPerHair, nEdges;
    static public float steph;

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
        //codeMatrixExample();
    }

    // Update is called once per frame
    void Update()
    {
        tickSimulation();
        updateClothVertexPositions();
    }

    void initData()
    {
        simulationSteps = 40;

        // hair nodes
        nHairs = 1;
        nNodesPerHair = 2;

        // global simulation variables
        nodeDistance = 0.5f;    // Initial Node distance apart.
        dPosition = 0.0004f;    // Euler method integration ratio for speed.
        dVelocity = 1.0f;       // Euler method integration ratio for acceleration.
        forceDecay = 0.0000f;
        velocityDecay = 0.999f;
        gravity = 0.1f;
        stiffness = 6.0f;
        maxTravelDistance = 5.0f;
        bendingStiffness = 0.1f;

        //cloth NodeList
        clothAttrsArray.pos = new Vector3[nHairs, nNodesPerHair];
        clothAttrsArray.vel = new Vector3[nHairs, nNodesPerHair];
        clothAttrsArray.acc = new Vector3[nHairs, nNodesPerHair];

        // fill data for each node
        int nodeIndex;
        double[] initial_pos = new double[nHairs * nNodesPerHair * 3]; // initial position holder 
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                // SoA CPU method Data struct
                clothAttrsArray.pos[i, j] = new Vector3(nodeDistance * (i - nHairs / 2), -nodeDistance * (j - nNodesPerHair / 2), 0.0f);
                clothAttrsArray.vel[i, j] = new Vector3(0.0f, 0.0f, 0.0f);
                clothAttrsArray.acc[i, j] = new Vector3(0.0f, 0.0f, 0.0f);
            }
        }
    }

    void initGeo()
    {
        // instantiate hair objects
        clothVertArray = new GameObject[nHairs * nNodesPerHair];
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                Vector3 location = clothAttrsArray.pos[i, j];
                var newitem = Instantiate(hairPrefab, location, Quaternion.identity);
                clothVertArray[i * nNodesPerHair + j] = newitem;
            }
        }
    }

    void updateClothVertexPositions()
    {
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                clothVertArray[i * nNodesPerHair + j].transform.position = clothAttrsArray.pos[i, j];
            }
        }
    }

    void tickSimulation()
    {
        explicitEulerSimulationMethod();
    }


    //////////////////////////// Cloth data run time update ////////////////////////////
    void explicitEulerSimulationMethod()
    {
        for (int i = 0; i < simulationSteps; i++)
        {
            explicitEuler_tickVel();
            explicitEuler_tickForce();
            explicitEuler_tickLight();
            explicitEuler_tickIntegration();
        }
    }

    // kernel Force in CPU implementation 
    void explicitEuler_tickForce()
    {
        // calculate F among the nodes
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair - 1; j++)
            {
                // spring force
                Vector3 currNext = clothAttrsArray.pos[i, j + 1] - clothAttrsArray.pos[i, j];
                float dX = Mathf.Clamp(nodeDistance - currNext.magnitude, -maxTravelDistance, maxTravelDistance);
                Vector3 springForce = -stiffness * currNext.normalized * dX;

                // bending force
                Vector3 bendingForce = new Vector3(0.0f, 0.0f, 0.0f);

                if (j != 0)
                {
                    Vector3 currPrev = clothAttrsArray.pos[i, j - 1] - clothAttrsArray.pos[i, j];
                    bendingForce = bendingStiffness * (currNext + currPrev);
                    clothAttrsArray.acc[i, j - 1] -= 0.5f * bendingForce;
                }

                clothAttrsArray.acc[i, j] += bendingForce;
                clothAttrsArray.acc[i, j] += springForce;
                clothAttrsArray.acc[i, j].y -= gravity;

                clothAttrsArray.acc[i, j + 1] -= springForce;
                clothAttrsArray.acc[i, j + 1] -= 0.5f * bendingForce;
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


    //////////////////////////// Slider Callbacks ////////////////////////////
    public void updateStiffness(float value)
    {
        stiffness = value;
    }

    public void updateGravity(float value)
    {
        gravity = value;
    }
}
