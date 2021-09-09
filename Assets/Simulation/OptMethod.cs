// Math Library
// https://numerics.mathdotnet.com/

// Math Plugin: https://github.com/GlitchEnzo/NuGetForUnity, following the steps to install MathNet:
// 1. add the NugetForUnity.unitypackage to your Unity Project.
// 2. Open the NuGet > Manage NuGet Packages Window.
// 3. Search the MathNet.Numerics and install.

// Mass-Spring System Solver
// https://github.com/alecjacobson/computer-graphics-mass-spring-systems?fbclid=IwAR1rxSAUgD4SzL6UOU4WBMKNZ8Rx5aKF-nNmpZ_Gx4w_fDHtckS-fKSXP6A
// 

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

public enum Selection { GPU, CPU, LocalGlobal };

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

// for Implicit Euler
public struct SpringConstraint
{
    // Spring constraints is for each spring, it is connected to two vertices, m_v1 and m_v2.
    public int m_v1, m_v2; 
    public float ks;
}

public struct PinConstraint
{
    // pin constraint is for each pinned vertex.
    public int m_v0;
}

public struct ClothSimulation
{
    public Matrix<double> rest_p; // 3m x 1
    public Matrix<double> prev_p; // 3m x 1
    public Matrix<double> curr_p; // 3m x 1

    public Matrix<double> inertia_y; // 3m x 1, y = 2 * curr_pos - prev_pos. 
    public Matrix<double> mass_matrix; // 3m x 3m, M with mass along the diagonal line. 
    public Matrix<double> external_force; // 3m x 1, F= ma, 
    public Matrix<double> gravity_force; // 3m x 1, F= -mg, external_force[3*i+1] is where gravity should apply

    public Matrix<double> A_matrix; // 3s x 3m
    public Matrix<double> l_matrix; // 3m x 3m
    public Matrix<double> j_matrix; // 3m x 3s

    public SpringConstraint[] springConstraint; // s, the number of springs
    public PinConstraint[] pinConstraint; // number of pinned vertices
}

public class OptMethod : MonoBehaviour
{
    // stuff can be seen from the UI
    public GameObject hairPrefab, colliderPrefab;
    public LampOptMethod theLamp;
    public Texture2D weightMap;
    public Selection simulationMode = Selection.LocalGlobal;

    // GPU Method Data Structure.
    static ClothNode[] clothNodesArray;
    static GameObject[] clothVertArray;

    // CPU Method Data Structure.
    ClothNodesAttrArray clothAttrsArray;

    // Implicit Method Data Structure.
    static ClothSimulation clothSim;

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

        // initialize empty data structures
        switch (simulationMode)
        {
            case Selection.GPU:
                clothNodesArray = new ClothNode[nHairs * nNodesPerHair];
                break;
            case Selection.CPU:
                //cloth NodeList
                clothAttrsArray.pos = new Vector3[nHairs, nNodesPerHair];
                clothAttrsArray.vel = new Vector3[nHairs, nNodesPerHair];
                clothAttrsArray.acc = new Vector3[nHairs, nNodesPerHair];
                break;
            case Selection.LocalGlobal:
                // implementation here
                nHairs = 10;
                nNodesPerHair = 10;
                nEdges = (nHairs - 1) * nNodesPerHair + (nNodesPerHair - 1) * nHairs;
                steph = 1.0f;
                nodeDistance = 3;
                stiffness = 5.0f;
                gravity = 0.5f;
                clothSim.springConstraint = new SpringConstraint[nEdges];
                clothSim.pinConstraint = new PinConstraint[2];
                break;
        }

        // fill data to Cloth Node Data structures
        int nodeIndex;
        double[] temp_p = new double[nHairs * nNodesPerHair * 3];
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {                
                switch (simulationMode)
                {
                    case Selection.GPU:
                        nodeIndex = i * nNodesPerHair + j;
                        clothNodesArray[nodeIndex].x = nodeDistance * (i - nHairs / 2);
                        clothNodesArray[nodeIndex].y = -nodeDistance * (j - nNodesPerHair / 2);
                        clothNodesArray[nodeIndex].z = 0.0f;
                        clothNodesArray[nodeIndex].vx = 0.0f;
                        clothNodesArray[nodeIndex].vy = 0.0f;
                        clothNodesArray[nodeIndex].vz = 0.0f;
                        clothNodesArray[nodeIndex].ax = 0;
                        clothNodesArray[nodeIndex].ay = 0;
                        clothNodesArray[nodeIndex].az = 0;

                        // sample the weight from the weight map
                        float u = i * 1.0f / (nHairs);
                        float v = 1.0f - j * 1.0f / (nNodesPerHair);
                        clothNodesArray[nodeIndex].mass = 1.0f + weightMap.GetPixelBilinear(u, v).grayscale;

                        break;
                    case Selection.CPU:
                        // SoA CPU method Data struct
                        clothAttrsArray.pos[i, j] = new Vector3(nodeDistance * (i - nHairs / 2), -nodeDistance * (j - nNodesPerHair / 2), 0.0f);
                        clothAttrsArray.vel[i, j] = new Vector3(0.0f, 0.0f, 0.0f);
                        clothAttrsArray.acc[i, j] = new Vector3(0.0f, 0.0f, 0.0f);
                        break;
                    case Selection.LocalGlobal:
                        // implementation here
                        nodeIndex = i * nNodesPerHair + j;

                        temp_p[3 * nodeIndex] =  nodeDistance * (i - nHairs / 2);
                        temp_p[3 * nodeIndex + 1] = -nodeDistance * (j - nNodesPerHair / 2);
                        temp_p[3 * nodeIndex + 2] = 0.0f;
                        break;
                }
            }
        }

        // filling Constraint Data Struct
        if (simulationMode == Selection.LocalGlobal)
        {
            // initialize clothMeshMatrix 
            clothSim.rest_p = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, temp_p);
            clothSim.prev_p = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, temp_p);
            clothSim.curr_p = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, temp_p);
            clothSim.inertia_y = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, temp_p);
            clothSim.mass_matrix = Matrix<double>.Build.DenseDiagonal(nHairs * nNodesPerHair * 3, nHairs * nNodesPerHair * 3, 1.0f);
            clothSim.external_force = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, 1, 0.0f);
            clothSim.gravity_force = buildGravityMatrix(gravity, nHairs * nNodesPerHair);
            clothSim.A_matrix = Matrix<double>.Build.Dense(nEdges * 3, nHairs * nNodesPerHair * 3, 0.0f);
            clothSim.l_matrix = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, nHairs * nNodesPerHair * 3, 0.0f);
            clothSim.j_matrix = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, nEdges * 3, 0.0f);

            // one spring
            clothSim.pinConstraint[0].m_v0 = 0;
            clothSim.pinConstraint[1].m_v0 = nNodesPerHair * (nHairs-1);

            int edgeIndex = 0;
            for (int i = 0; i < nHairs ; i++)
            {
                for (int j = 0; j < nNodesPerHair; j++)
                {
                    int vertexIndex = i * nNodesPerHair + j;

                    // horizontal constraint
                    if (i + 1 < nHairs) {
                        clothSim.springConstraint[edgeIndex].m_v1 = vertexIndex;
                        clothSim.springConstraint[edgeIndex].m_v2 = vertexIndex + nNodesPerHair;
                        clothSim.springConstraint[0].ks = stiffness;
                        //Debug.Log(">>> " + edgeIndex + ": " + clothSim.springConstraint[edgeIndex].m_v1 + ", " + clothSim.springConstraint[edgeIndex].m_v2);
                        edgeIndex += 1;
                    }

                    if (j + 1 < nNodesPerHair)
                    {
                        clothSim.springConstraint[edgeIndex].m_v1 = vertexIndex;
                        clothSim.springConstraint[edgeIndex].m_v2 = vertexIndex + 1;
                        clothSim.springConstraint[0].ks = stiffness;
                        //Debug.Log(">>> " + edgeIndex + ": " + clothSim.springConstraint[edgeIndex].m_v1 + ", " + clothSim.springConstraint[edgeIndex].m_v2);
                        edgeIndex += 1;
                    }
                }
            }

            //clothSim.springConstraint[0].m_v1 = 0;
            //clothSim.springConstraint[0].m_v2 = 1;
            //clothSim.springConstraint[0].ks = stiffness;

            //// second spring
            //clothSim.springConstraint[1].m_v1 = 1;
            //clothSim.springConstraint[1].m_v2 = 2;
            //clothSim.springConstraint[1].ks = stiffness;

            //// third spring
            //clothSim.springConstraint[2].m_v1 = 2;
            //clothSim.springConstraint[2].m_v2 = 3;
            //clothSim.springConstraint[2].ks = stiffness;

            //clothSim.springConstraint[3].m_v1 = 3;
            //clothSim.springConstraint[3].m_v2 = 4;
            //clothSim.springConstraint[3].ks = stiffness;
        }
    }

    void initGeo()
    {
        // instantiate hair objects
        clothVertArray = new GameObject[nHairs * nNodesPerHair];
        switch (simulationMode)
        {
            case Selection.GPU:
                for (int i = 0; i < clothVertArray.Length; i++)
                {
                    Vector3 location = new Vector3(clothNodesArray[i].x, clothNodesArray[i].y, clothNodesArray[i].z);
                    var newitem = Instantiate(hairPrefab, location, Quaternion.identity);
                    clothVertArray[i] = newitem;
                }
                break;
            case Selection.CPU:
                for (int i = 0; i < nHairs; i++)
                {
                    for (int j = 0; j < nNodesPerHair; j++)
                    {
                        Vector3 location = clothAttrsArray.pos[i, j];
                        var newitem = Instantiate(hairPrefab, location, Quaternion.identity);
                        clothVertArray[i * nNodesPerHair + j] = newitem;
                    }
                }
                break;
            case Selection.LocalGlobal:
                // your implementation here
                for (int i = 0; i < nHairs * nNodesPerHair; i++)
                {
                    Vector3 location = new Vector3(Convert.ToSingle(clothSim.curr_p.At(3 * i, 0)),
                                                   Convert.ToSingle(clothSim.curr_p.At(3 * i + 1, 0)),
                                                   Convert.ToSingle(clothSim.curr_p.At(3 * i + 2, 0)));
                    GameObject newitem = Instantiate(hairPrefab, location, Quaternion.identity);
                    newitem.transform.localScale = new Vector3(1, 1, 1);
                    clothVertArray[i] = newitem;
                }
                break;
        }

    }

    void updateClothVertexPositions()
    {
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                switch (simulationMode)
                {
                    case Selection.GPU:
                        // implementation
                        break;
                    case Selection.CPU:
                        clothVertArray[i * nNodesPerHair + j].transform.position = clothAttrsArray.pos[i, j];
                        break;
                    case Selection.LocalGlobal:
                        // implementation
                        int nodeIndex = i * nNodesPerHair + j;
                        clothVertArray[nodeIndex].transform.position = new Vector3(Convert.ToSingle(clothSim.curr_p.At(3 * nodeIndex, 0)),
                                                                                   Convert.ToSingle(clothSim.curr_p.At(3 * nodeIndex + 1, 0)),
                                                                                   Convert.ToSingle(clothSim.curr_p.At(3 * nodeIndex + 2, 0)));
                        break;
                }                
            }
        }
    }

    void tickSimulation()
    {
        switch (simulationMode)
        {
            case Selection.GPU:
                GPUSimulationMethod();
                break;
            case Selection.CPU:
                explicitEulerSimulationMethod();
                // p [matrix of vector3 # of nodes I have] = p+v*h
                // v = v+ a*m
                // for every node. 
                break;
            case Selection.LocalGlobal:
                localGlobalSimulationMethod();
                // finish implementation? understand what is going on here
                // Expectation: continue this research. 
                // by early fall, be able 
                break;
            default:
                Debug.Log("Unknown Simulation Method!");
                break;
        }
    }

    void localGlobalSimulationMethod()
    {
        // Step1 CalculateInteria component y = 2 * P_current - P_previous
        LocalGlobal_updateInertiaY();      

        // Step2 build External Force Vector
        LocalGlobal_updateExternalForce();

        // Calculate L matrix, LLaplacian generate smooth transition from A to B. Investigate this.
        LocalGlobal_updateAandLMatrix();

        // Calculate J matrix
        LocalGlobal_updateJMatrix();
        clothSim.j_matrix = stiffness * clothSim.A_matrix.Transpose();

        // Calculate Pinned constraint CT * C
        Matrix<double> CTC = LocalGlobal_evaluate_CTC(10000000);

        // evaluate Q = M + h^2 * ( L + CTC )        
        Matrix<double> Q = clothSim.mass_matrix + steph * steph * (clothSim.l_matrix + CTC) ;

        // evaluate d
        Matrix<double> d = LocalGlobal_evaluate_d();

        // calculate b = M * y + (f + J * d + CTC * clothSim.rest_p) * h^2
        Matrix<double> b = clothSim.mass_matrix * clothSim.inertia_y + (clothSim.external_force + clothSim.j_matrix * d + CTC * clothSim.rest_p) * steph * steph;

        // solve equation Q p = b, Q = LT * L
        Matrix<double> next_p = Q.Solve(b);
        next_p.SetSubMatrix(0, 3, 0, 1, clothSim.rest_p.SubMatrix(0, 3, 0, 1));

        // upgrade current, previous position. next_p.SetSubMatrix(0, 3, 0, 1, clothSim.rest_p.SubMatrix(0, 3, 0, 1));
        clothSim.prev_p = clothSim.curr_p;
        clothSim.curr_p = next_p;

        // Debug.Log(clothSim.curr_p);

        // calculating damping

        return;
    }

    void LocalGlobal_updateInertiaY()
    {
        //y = 2 * curr_pos - prev_pos.
        clothSim.inertia_y = clothSim.curr_p + (clothSim.curr_p - clothSim.prev_p);
    }

    void LocalGlobal_updateExternalForce()
    {
        clothSim.external_force = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, 1, 0.0f);
        clothSim.gravity_force = buildGravityMatrix(gravity, nHairs * nNodesPerHair);
        clothSim.external_force = clothSim.external_force + clothSim.gravity_force;
        clothSim.external_force = clothSim.mass_matrix * clothSim.external_force;
        //Debug.Log(clothSim.external_force);
    }

    void LocalGlobal_updateLaplacian()
    {
        // This implementation is wrong.
        // L is 3m * 3m size
        double sk = Convert.ToDouble(stiffness);
        Tuple<int, int, double>[] L = new Tuple<int, int, double>[12 * clothSim.springConstraint.Length];
        for (int i = 0; i < clothSim.springConstraint.Length; i++)
        {
            SpringConstraint s = clothSim.springConstraint[i];
            L[12 * i + 0] = Tuple.Create(3 * s.m_v1 + 0, 3 * s.m_v1 + 0, sk);
            L[12 * i + 1] = Tuple.Create(3 * s.m_v1 + 1, 3 * s.m_v1 + 1, sk);
            L[12 * i + 2] = Tuple.Create(3 * s.m_v1 + 2, 3 * s.m_v1 + 2, sk);

            L[12 * i + 3] = Tuple.Create(3 * s.m_v1 + 0, 3 * s.m_v2 + 0, -sk);
            L[12 * i + 4] = Tuple.Create(3 * s.m_v1 + 1, 3 * s.m_v2 + 1, -sk);
            L[12 * i + 5] = Tuple.Create(3 * s.m_v1 + 2, 3 * s.m_v2 + 2, -sk);

            L[12 * i + 6] = Tuple.Create(3 * s.m_v2 + 0, 3 * s.m_v1 + 0, -sk);
            L[12 * i + 7] = Tuple.Create(3 * s.m_v2 + 1, 3 * s.m_v1 + 1, -sk);
            L[12 * i + 8] = Tuple.Create(3 * s.m_v2 + 2, 3 * s.m_v1 + 2, -sk);

            L[12 * i + 9] = Tuple.Create(3 * s.m_v2 + 0, 3 * s.m_v2 + 0, sk);
            L[12 * i + 10] = Tuple.Create(3 * s.m_v2 + 1, 3 * s.m_v2 + 1, sk);
            L[12 * i + 11] = Tuple.Create(3 * s.m_v2 + 2, 3 * s.m_v2 + 2, sk);
        }

        clothSim.l_matrix = Matrix<double>.Build.DenseOfIndexed(nHairs * nNodesPerHair * 3, nHairs * nNodesPerHair * 3, L);
    }

    void LocalGlobal_updateAandLMatrix()
    {
        // L is 3m * 3m size
        double sk = Convert.ToDouble(stiffness);
        Tuple<int, int, double>[] A = new Tuple<int, int, double>[6 * clothSim.springConstraint.Length];
        for (int index = 0; index < clothSim.springConstraint.Length; index++)
        {
            SpringConstraint s = clothSim.springConstraint[index];
            A[6 * index + 0] = Tuple.Create(3 * index + 0, 3 * s.m_v1 + 0, 1.0);
            A[6 * index + 1] = Tuple.Create(3 * index + 1, 3 * s.m_v1 + 1, 1.0);
            A[6 * index + 2] = Tuple.Create(3 * index + 2, 3 * s.m_v1 + 2, 1.0);
            A[6 * index + 3] = Tuple.Create(3 * index + 0, 3 * s.m_v2 + 0, -1.0);
            A[6 * index + 4] = Tuple.Create(3 * index + 1, 3 * s.m_v2 + 1, -1.0);
            A[6 * index + 5] = Tuple.Create(3 * index + 2, 3 * s.m_v2 + 2, -1.0);
        }
        //Debug.Log(A.Length);
        //Debug.Log(clothSim.springConstraint.Length * 3);
        //Debug.Log(nHairs * nNodesPerHair * 3);
        //for (int i = 0; i < clothSim.springConstraint.Length; i++)
        //{
        //    Debug.Log(A[6 * i + 0] + ", " + A[6 * i + 1] + ", " + A[6 * i + 2] + ", ");
        //    Debug.Log(A[6 * i + 3] + ", " + A[6 * i + 4] + ", " + A[6 * i + 5] + ", ");
        //}
        clothSim.A_matrix = Matrix<double>.Build.DenseOfIndexed(clothSim.springConstraint.Length * 3, nHairs * nNodesPerHair * 3, A);
        Matrix<double> L = sk * clothSim.A_matrix.TransposeThisAndMultiply(clothSim.A_matrix);
        clothSim.l_matrix = L;
    }

    void LocalGlobal_updateJMatrix()
    {
        // J is 3m * 3s size
        double sk = Convert.ToDouble(stiffness);
        Tuple<int, int, double>[] J = new Tuple<int, int, double>[6 * clothSim.springConstraint.Length];
        for (int i = 0; i < clothSim.springConstraint.Length; i++)
        {
            SpringConstraint s = clothSim.springConstraint[i];
            J[6 * i + 0] = Tuple.Create(3 * s.m_v1 + 0, 3 * i + 0, sk);
            J[6 * i + 1] = Tuple.Create(3 * s.m_v1 + 1, 3 * i + 1, sk);
            J[6 * i + 2] = Tuple.Create(3 * s.m_v1 + 2, 3 * i + 2, sk);

            J[6 * i + 3] = Tuple.Create(3 * s.m_v2 + 0, 3 * i + 0, -sk);
            J[6 * i + 4] = Tuple.Create(3 * s.m_v2 + 1, 3 * i + 1, -sk);
            J[6 * i + 5] = Tuple.Create(3 * s.m_v2 + 2, 3 * i + 2, -sk);
        }
        clothSim.j_matrix = Matrix<double>.Build.DenseOfIndexed(nHairs * nNodesPerHair * 3, nEdges * 3, J);        
    }

    Matrix<double> LocalGlobal_evaluate_CTC(double penalty)
    {
        // CTC is 3m * 3m size
        Tuple<int, int, double>[] CTC = new Tuple<int, int, double>[3 * clothSim.pinConstraint.Length];
        for (int i = 0; i < clothSim.pinConstraint.Length; i++)
        {
            PinConstraint pin = clothSim.pinConstraint[i];
            CTC[3 * i + 0] = Tuple.Create(3 * pin.m_v0 + 0, 3 * pin.m_v0 + 0, penalty);
            CTC[3 * i + 1] = Tuple.Create(3 * pin.m_v0 + 1, 3 * pin.m_v0 + 1, penalty);
            CTC[3 * i + 2] = Tuple.Create(3 * pin.m_v0 + 2, 3 * pin.m_v0 + 2, penalty);
        }
        return Matrix<double>.Build.DenseOfIndexed(nHairs * nNodesPerHair * 3, nHairs * nNodesPerHair * 3, CTC);
    }

    Matrix<double> LocalGlobal_evaluate_d()
    {
        Matrix<double> d = Matrix<double>.Build.Dense(nEdges * 3, 1, 0.0f);
        for (int i = 0; i < nEdges; i++)
        {
            SpringConstraint s = clothSim.springConstraint[i];
            Matrix<double> p1 = clothSim.curr_p.SubMatrix(3*s.m_v1, 3, 0, 1);
            Matrix<double> p2 = clothSim.curr_p.SubMatrix(3*s.m_v2, 3, 0, 1);
            Matrix<double> di = (p1 - p2).NormalizeColumns(2.0) * nodeDistance;
            d.SetSubMatrix(3*i, 3, 0, 1, di);
        }
        return d;
    }

    Matrix<double> buildGravityMatrix(float gravity, int mSize)
    {
        // gravity is a 3m x 3m size, which only affects Y direction
        double[] temp_g = new double[3 * mSize];
        for (int i = 0; i < mSize; i++)
        {
            temp_g[3 * i + 1] = -gravity;
        }
        return Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, temp_g);
    }

    void explicitEulerSimulationMethod()
    {
        for (int i = 0; i < simulationSteps; i++)
        {            
            explicitEuler_tickVel();
            explicitEuler_tickForce();
            explicitEuler_tickLight();
            explicitEuler_tickIntegration();
            // do it in dynamic way, so that we don't have to through all iterations to reach destinations.
            // difference is small quit iteration. adjustment is small, quit. 
            // skip nodes. 
            // check if it should continue. 
            // define what is close?
            // do it linear 
            // integration may oscillate. Detect error/ oscillate, in the code, 
            // analysis, understand the relation to localGlobalSimulationMethod(); 
        }
    }

    void codeMatrixExample()
    {
        // basic operation
        Matrix<double> A = DenseMatrix.OfArray(new double[,] {  {1,1},
                                                                {1,2} });
        Matrix<double> B = DenseMatrix.OfArray(new double[,] {  {-1,-1},
                                                                {3, 4} });
        Debug.Log(A + B);

        // build Matrix from Vector
        double[] v = { 1, 2, 3 ,4};
        Matrix<double> C = Matrix<double>.Build.DenseOfColumnMajor(2, 2, v);
        Matrix<double> D = Matrix<double>.Build.DenseOfColumnMajor(4, 1, v);
        Matrix<double> E = Matrix<double>.Build.DenseOfColumnMajor(1, 4, v);
        Debug.Log(C);
        Debug.Log(D);
        Debug.Log(E);

        // build from index
        Tuple<int, int, double>[] t = { Tuple.Create(0, 0, 2.0), Tuple.Create(0, 1, -3.0) };
        Matrix<double> F = Matrix<double>.Build.DenseOfIndexed(3, 4, t);
        Debug.Log(F);

        // Matrix Operation
        Matrix<double> DE = D * E;
        Matrix<double> FD = F * D;
        Debug.Log(DE);
        Debug.Log(FD);

        // Modify the elements of a Matrix
        Matrix<double> SUB = Matrix<double>.Build.DenseOfArray(new double[,] { { 1983 } });
        DE.SetSubMatrix(2, 1, 2, 1, SUB);
        Debug.Log(DE);

        // solve equation Hx=b
        var H = Matrix<double>.Build.DenseOfArray(new double[,] {   { 3, 2, -1 },
                                                                    { 2, -2, 4 },
                                                                    { -1, 0.5, -1 } });
        var b = Vector<double>.Build.Dense(new double[] { 1, -2, 0 });
        var b2 = Matrix<double>.Build.DenseOfArray(new double[,] {  { 1 },
                                                                    { -2 },
                                                                    { 0 } });
        var x = H.Solve(b2);
        Debug.Log(x);
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
