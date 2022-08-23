// Math Library documentation
// https://numerics.mathdotnet.com/
// https://numerics.mathdotnet.com/api/MathNet.Numerics.LinearAlgebra.Complex/Matrix.htm#SetSubMatrix

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

// Implicit Euler Spring Constraint
public struct SpringConstraint
{
    // Spring constraints is for each spring, it is connected to two vertices, m_v1 and m_v2.
    public int m_v1, m_v2; 
    public float ks;
}

// Implicit Euler Pin Constraint
public struct PinConstraint
{
    // pin constraint is for each pinned vertex.
    public int m_v0;
}

public struct ColliderNodes
{
    public float x, y, z, r;
};

public struct ClothSimulationImplicit
{
    // m is size of the nodes, s is the size of the spring. 
    public Matrix<double> rest_p; // 3m x 1, flattend (x, y, z) collumn matrix. 
    public Matrix<double> prev_p; // 3m x 1
    public Matrix<double> curr_p; // 3m x 1

    public Matrix<double> inertia_y; // 3m x 1, y = 2 * curr_pos - prev_pos. 
    public Matrix<double> mass_matrix; // 3m x 3m, M with mass along the diagonal line. 
    
    public Matrix<double> gravity_unit; // 3m x 1, gravity_unit[3*1+1] is to normalized constant 1;
    public Matrix<double> gravity_force; // 3m x 1, gravity_force = g * gravity_unit, this matrix is updated in run time.
    public Matrix<double> light_force; // light cast on the nodes, generating a force matrix. 
    public Matrix<double> external_force; // 3m x 1, F = gravity + light.

    public Matrix<double> A_matrix; // 3s x 3m, A defines the relations between Spring and Nodes. A_ij = 1, -1 or 0;
    public Matrix<double> AT_matrix; // 3m X 3s, AT = A.transpose();
    public Matrix<double> ATA_matrix; // 3m x 3m, ATA = AT * A, ATA should be pre-computed to update L matrix;
    public Matrix<double> L_matrix; // 3m x 3m, L is the stiffness weighted Laplacian among the Nodes. L = k * ATA;
    public Matrix<double> J_matrix; // 3m x 3s, J = k * AT;
    public Matrix<double> CTC_matrix; // 3m x 3m, CTC = C.transpose() * C. C this the contraint matrix;

    public SpringConstraint[] springConstraint; // size = the number of springs
    public PinConstraint[] pinConstraint; // size = number of pinned vertices
}

public class OptMethod : MonoBehaviour
{
    public GlobalData globalData;

    // stuff can be seen from the UI
    public GameObject hairPrefab, sphereColliderPrefab, planeColliderPrefab;
    public Lamp theLamp;
    public Texture2D weightMap;
    public int LightIntensityMutiplier = 10;
    public int resolution = 10;
    public SimulationScenario simulationScenario;
    public ClothMeshGenerator clothMesh;
    public float clothDebugNodeSize;
    public FlagController theFlag;

    public float nodeDistance;
    public float gravity;
    public float stiffness;

    // Cloth vertex Data Structure.
    static GameObject[] clothVertArray;

    // Implicit Method Data Structure.
    static ClothSimulationImplicit ClothSimImp;

    // hair nodes, hair is a vertical line of the cloth
    static public int nHairs, nNodesPerHair, nEdges;
    static public float steph;

    // sphere colliders
    public float collisionPushAwayDistance;
    static public int nColliders;
    static public float colliderRadius;
    static ColliderNodes[] colliderNodeArrays;
    static public GameObject[] colliderGeos;

    // plane collider
    static public float planeColliderPosY = 0;

    public float lightForce;
    public float lightAngle;
    public Vector3 headPos;
    public Vector3 headLookAtPos;

    // pin constraint penalty amount
    private float PIN_PENALTY = 10000000;

    // Start is called before the first frame update
    void Start()
    {
        Debug.Assert(globalData);
        Debug.Assert(hairPrefab);
        Debug.Assert(sphereColliderPrefab);
        Debug.Assert(planeColliderPrefab);
        Debug.Assert(theLamp);
        Debug.Assert(weightMap);
        Debug.Assert(clothMesh);
        Debug.Assert(theFlag);
        Debug.Assert(resolution > 0);
        Debug.Assert(nodeDistance > 0);

        initData();
        initClothVertices();
        initClothMesh();
        initFlagController();
    }

    // Update is called once per frame
    void Update()
    {
        tickSimulation();
        updateClothVertices();
        updateClothMesh();
    }

    //////////////////////////// Cloth Initialization ////////////////////////////
    void initClothVertices()
    {
        // instantiate hair objects
        clothVertArray = new GameObject[nHairs * nNodesPerHair];
        for (int i = 0; i < nHairs * nNodesPerHair; i++)
        {
            Vector3 location = getVector3FromColumnMatrix(ClothSimImp.curr_p, i);
            GameObject newitem = Instantiate(hairPrefab, location, Quaternion.identity);
            newitem.transform.localScale = new Vector3(1, 1, 1) * clothDebugNodeSize;
            newitem.GetComponent<NodeController>().setNodeIndex(i);
            clothVertArray[i] = newitem;
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
        clothMesh.initMeshImplicit(ClothSimImp.curr_p, nHairs, nNodesPerHair);
    }

    void updateClothMesh()
    {
        clothMesh.updateMeshImplict(ClothSimImp.curr_p);
    }

    void initFlagController()
    {
        theFlag.gameObject.SetActive(false);
    }

    void initData()
    {
        clothDebugNodeSize = globalData.clothDebugNodeSize;

        // hair nodes
        nHairs = globalData.resolution;
        nNodesPerHair = globalData.resolution;
        nEdges = (nHairs - 1) * nNodesPerHair + (nNodesPerHair - 1) * nHairs;

        // global simulation variables
        nodeDistance = globalData.nodeDistance;
        gravity = globalData.gravity;
        stiffness = globalData.stiffness;
        steph = 1.0f;

        // light control
        lightForce = globalData.lightForce;
        lightAngle = globalData.lightAngle;

        headPos = theLamp.head.transform.position;
        headLookAtPos = theLamp.headLookAt;

        // empty constraints
        ClothSimImp.springConstraint = new SpringConstraint[nEdges];
        ClothSimImp.pinConstraint = new PinConstraint[2];

        // fill data for each node
        int nodeIndex;
        double[] initial_pos = new double[nHairs * nNodesPerHair * 3]; // initial position holder 
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {                
                // implementation here
                nodeIndex = i * nNodesPerHair + j;
                if (simulationScenario == SimulationScenario.Static)
                {
                    initial_pos[3 * nodeIndex] = nodeDistance * (i - nHairs / 2);
                    initial_pos[3 * nodeIndex + 1] = -nodeDistance * (j - nNodesPerHair / 2);
                    initial_pos[3 * nodeIndex + 2] = 0.0f;
                }

                if (simulationScenario == SimulationScenario.FreeFall)
                {
                    initial_pos[3 * nodeIndex] = nodeDistance * (i - nHairs / 2);
                    initial_pos[3 * nodeIndex + 1] = 0.0f;
                    initial_pos[3 * nodeIndex + 2] = -nodeDistance * (j - nNodesPerHair / 2);
                }

                if (simulationScenario == SimulationScenario.Flag)
                {
                    initial_pos[3 * nodeIndex] = nodeDistance * (i - nHairs / 2);
                    initial_pos[3 * nodeIndex + 1] = -nodeDistance * (j - nNodesPerHair / 2);
                    initial_pos[3 * nodeIndex + 2] = 0.0f;
                }

            }
        }

        // initialize clothMeshMatrix 
        ClothSimImp.rest_p = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, initial_pos);
        ClothSimImp.prev_p = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, initial_pos);
        ClothSimImp.curr_p = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, initial_pos);
        ClothSimImp.inertia_y = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, initial_pos);
        ClothSimImp.mass_matrix = Matrix<double>.Build.DenseDiagonal(nHairs * nNodesPerHair * 3, nHairs * nNodesPerHair * 3, 1.0f);
        ClothSimImp.external_force = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, 1, 0.0f);
        ClothSimImp.light_force = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, 1, 0.0f) ;
        ClothSimImp.gravity_unit = LocalGlobal_compute_MatrixG(1, nHairs * nNodesPerHair);
        ClothSimImp.A_matrix = Matrix<double>.Build.Dense(nEdges * 3, nHairs * nNodesPerHair * 3, 0.0f);
        ClothSimImp.L_matrix = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, nHairs * nNodesPerHair * 3, 0.0f);
        ClothSimImp.J_matrix = Matrix<double>.Build.Dense(nHairs * nNodesPerHair * 3, nEdges * 3, 0.0f);

        // set the pinned nodes
        ClothSimImp.pinConstraint[0].m_v0 = 0;
        ClothSimImp.pinConstraint[1].m_v0 = nNodesPerHair * (nHairs-1);

        // set up mass-srping system
        // 1. every spring is connected by two nodes, v1 and v2. 
        // 2. every spring has a stiffness coefficient. 
        // 3. every spring is indexed.
        // 
        //             idx:0        idx:n      idx:2n
        //       (V1) ------- (V2) ------- () ------- 
        //        |            |
        //        |  idx:1     |  idx:n+1
        //        |
        //       (V2')------- ( ) ---
        //        |            |
        //        |  idx:2     |
        //        |
        //       ( ) ---

        int edgeIndex = 0;
        for (int i = 0; i < nHairs ; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                int vertexIndex = i * nNodesPerHair + j;

                // Horizontal spring constraints.
                if (i + 1 < nHairs) {
                    ClothSimImp.springConstraint[edgeIndex].m_v1 = vertexIndex;
                    ClothSimImp.springConstraint[edgeIndex].m_v2 = vertexIndex + nNodesPerHair;
                    ClothSimImp.springConstraint[edgeIndex].ks = stiffness;
                    edgeIndex += 1;
                }

                // vertical spring constraint
                if (j + 1 < nNodesPerHair)
                {
                    ClothSimImp.springConstraint[edgeIndex].m_v1 = vertexIndex;
                    ClothSimImp.springConstraint[edgeIndex].m_v2 = vertexIndex + 1;
                    ClothSimImp.springConstraint[edgeIndex].ks = stiffness;
                    edgeIndex += 1;
                }
            }
        }

        // pre-compute A, A transpose, and then AT*A matrix
        LocalGlobal_compute_matrixA();
        LocalGlobal_compute_matrixAT();
        LocalGlobal_compute_matrixATA();

        // calculate L matrix, L =  k * ATA;
        LocalGlobal_update_matrixL();

        // calculate CTC matrix, CTC is used to set the pinned vertex;
        LocalGlobal_compute_MatrixCTC(PIN_PENALTY);

        // sphere collider data
        nColliders = 3;
        colliderRadius = 4.0f;
        colliderNodeArrays = new ColliderNodes[nColliders];
        for (int i = 0; i < nColliders; i++)
        {
            colliderNodeArrays[i].x = colliderRadius * (i - nColliders / 2); ;
            colliderNodeArrays[i].y = -20.0f;
            colliderNodeArrays[i].z = 0.0f;
            colliderNodeArrays[i].r = colliderRadius;
        }

        // plane collider data
        planeColliderPosY = -1 * resolution * nodeDistance / 2.0f;
    }

    //////////////////////////// Cloth Simulation Loop ////////////////////////////
    void tickSimulation()
    {
        implicitEulerSimulationMethod();
    }

    void implicitEulerSimulationMethod()
    {
        // CalculateInteria component y = 2 * P_current - P_previous
        LocalGlobal_update_inertiaY();      

        // Build External Force Vector
        LocalGlobal_update_externalForce();

        // Calculate L matrix, L = k * ATA; 
        LocalGlobal_update_matrixL();

        // Calculate J matrix, J = K * AT;
        LocalGlobal_update_matrixJ();

        // evaluate Q = M + h^2 * ( L + CTC )        
        Matrix<double> Q = ClothSimImp.mass_matrix + steph * steph * (ClothSimImp.L_matrix + ClothSimImp.CTC_matrix) ;

        // evaluate d
        Matrix<double> d = LocalGlobal_evaluate_d();

        // calculate b = M * y + (f + J * d + CTC * clothSim.rest_p) * h^2
        Matrix<double> b = ClothSimImp.mass_matrix * ClothSimImp.inertia_y + (ClothSimImp.external_force + ClothSimImp.J_matrix * d + ClothSimImp.CTC_matrix * ClothSimImp.rest_p) * steph * steph;

        // Find the next position by solving equation Q * p = b. To Do: prefactor Q = LT * L.
        Matrix<double> next_p = Q.Solve(b);

        // upgrade current, previous position.
        ClothSimImp.prev_p = ClothSimImp.curr_p;
        ClothSimImp.curr_p = next_p;

        // collision
        LocalGlobal_update_position_from_collision();

        // calculating damping?
    }

    //////////////////////////// Cloth data pre-computation ////////////////////////////
    void LocalGlobal_compute_matrixA()
    {
        Tuple<int, int, double>[] A = new Tuple<int, int, double>[6 * ClothSimImp.springConstraint.Length];
        for (int index = 0; index < ClothSimImp.springConstraint.Length; index++)
        {
            SpringConstraint s = ClothSimImp.springConstraint[index];
            A[6 * index + 0] = Tuple.Create(3 * index + 0, 3 * s.m_v1 + 0, 1.0);
            A[6 * index + 1] = Tuple.Create(3 * index + 1, 3 * s.m_v1 + 1, 1.0);
            A[6 * index + 2] = Tuple.Create(3 * index + 2, 3 * s.m_v1 + 2, 1.0);
            A[6 * index + 3] = Tuple.Create(3 * index + 0, 3 * s.m_v2 + 0, -1.0);
            A[6 * index + 4] = Tuple.Create(3 * index + 1, 3 * s.m_v2 + 1, -1.0);
            A[6 * index + 5] = Tuple.Create(3 * index + 2, 3 * s.m_v2 + 2, -1.0);
        }
        ClothSimImp.A_matrix = Matrix<double>.Build.DenseOfIndexed(ClothSimImp.springConstraint.Length * 3, nHairs * nNodesPerHair * 3, A);
    }

    void LocalGlobal_compute_MatrixCTC(double penalty)
    {
        // CTC is 3m * 3m size
        Tuple<int, int, double>[] CTC = new Tuple<int, int, double>[3 * ClothSimImp.pinConstraint.Length];
        for (int i = 0; i < ClothSimImp.pinConstraint.Length; i++)
        {
            PinConstraint pin = ClothSimImp.pinConstraint[i];
            CTC[3 * i + 0] = Tuple.Create(3 * pin.m_v0 + 0, 3 * pin.m_v0 + 0, penalty);
            CTC[3 * i + 1] = Tuple.Create(3 * pin.m_v0 + 1, 3 * pin.m_v0 + 1, penalty);
            CTC[3 * i + 2] = Tuple.Create(3 * pin.m_v0 + 2, 3 * pin.m_v0 + 2, penalty);
        }
        ClothSimImp.CTC_matrix = Matrix<double>.Build.DenseOfIndexed(nHairs * nNodesPerHair * 3, nHairs * nNodesPerHair * 3, CTC);
    }

    void setCTCMatrixAtIndex(int m_v0, float penalty)
    {
        Matrix<double> subPinConstraintMatrix = Matrix<double>.Build.DenseOfArray(new double[,] {  { penalty, 0, 0 },
                                                                                                 { 0, penalty, 0 },
                                                                                                 { 0, 0, penalty }  });
        ClothSimImp.CTC_matrix.SetSubMatrix(3 * m_v0, 3, 3 * m_v0, 3, subPinConstraintMatrix); // void SetSubMatrix(int rowIndex, int rowCount, int columnIndex, int columnCount, Matrix<T> subMatrix)
    }

    void setRestPositionAtIndex(int index)
    {
        Vector3 currentPosAtIndex = getVector3FromColumnMatrix(ClothSimImp.curr_p, index);
        setColumnMatrixFromVector3(ClothSimImp.rest_p, index, currentPosAtIndex);
    }

    public void setHairNodesPinStatusAtIndex(int nodeIndex, int status)
    {
        setCTCMatrixAtIndex(nodeIndex, status * PIN_PENALTY);
        setRestPositionAtIndex(nodeIndex);
    }

    Matrix<double> LocalGlobal_compute_MatrixG(float gravity, int mSize)
    {
        // gravity matrix is a 3m x 3m size, which only affects in the Y direction
        double[] temp_g = new double[3 * mSize];
        for (int i = 0; i < mSize; i++)
        {
            temp_g[3 * i + 1] = -gravity;
        }
        return Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, temp_g);
    }

    void LocalGlobal_compute_matrixAT()
    { 
        ClothSimImp.AT_matrix = ClothSimImp.A_matrix.Transpose();
    }

    void LocalGlobal_compute_matrixATA()
    {
        ClothSimImp.ATA_matrix = ClothSimImp.A_matrix.TransposeThisAndMultiply(ClothSimImp.A_matrix);
    }


    //////////////////////////// Cloth data run time update ////////////////////////////
    void updateClothVertices()
    {
        for (int i = 0; i < nHairs; i++)
        {
            for (int j = 0; j < nNodesPerHair; j++)
            {
                // implementation
                int nodeIndex = i * nNodesPerHair + j;
                clothVertArray[nodeIndex].transform.position = getVector3FromColumnMatrix(ClothSimImp.curr_p, nodeIndex);
            }
        }
    }

    void LocalGlobal_update_position_from_collision()
    {
        // double[] new_pos = new double[nHairs * nNodesPerHair * 3];      
        for (int j = 0; j < nColliders; j++)
        {
            for (int i = 0; i < nHairs * nNodesPerHair; i++)
            {
                Vector3 curr_node_pos = getVector3FromColumnMatrix(ClothSimImp.curr_p, i);
                Vector3 node_to_collider = curr_node_pos - colliderGeos[j].transform.position;
                colliderRadius = colliderGeos[j].transform.localScale.x / 2.0f;
                float totalPushAwayDistance = colliderRadius + collisionPushAwayDistance;

                // sphere collider
                if (node_to_collider.magnitude < totalPushAwayDistance) {
                    Vector3 p = colliderGeos[j].transform.position + node_to_collider.normalized * totalPushAwayDistance;
                    setColumnMatrixFromVector3(ClothSimImp.curr_p, i, p);
                    setColumnMatrixFromVector3(ClothSimImp.prev_p, i, p);
                    setColumnMatrixFromVector3(ClothSimImp.rest_p, i, p);
                }

                // plane collider
                if (curr_node_pos.y < planeColliderPosY)
                {
                    Vector3 p = new Vector3(curr_node_pos.x, planeColliderPosY, curr_node_pos.z);
                    setColumnMatrixFromVector3(ClothSimImp.curr_p, i, p);
                    setColumnMatrixFromVector3(ClothSimImp.prev_p, i, p);
                    setColumnMatrixFromVector3(ClothSimImp.rest_p, i, p);
                }
            }
        }
    return;
    }

    void LocalGlobal_update_inertiaY()
    {
        // y = 2 * curr_pos - prev_pos.
        ClothSimImp.inertia_y = ClothSimImp.curr_p + (ClothSimImp.curr_p - ClothSimImp.prev_p);
    }

    void LocalGlobal_update_lightForce() {
        
        double[] light_force_array = new double[nHairs * nNodesPerHair * 3];

        // lookAtDir
        Vector3 head_look_at_dir = (headLookAtPos - headPos).normalized;

        for (int i = 0; i < nHairs * nNodesPerHair; i++) {
            Vector3 curr_node_pos = getVector3FromColumnMatrix(ClothSimImp.curr_p, i);
            Vector3 head_to_node_dir = (curr_node_pos - headPos).normalized;
            float dot_value = Vector3.Dot(head_look_at_dir, head_to_node_dir);

            if (dot_value > Math.Cos(lightAngle * 0.0174533f * 0.5f)) {
                light_force_array[3 * i + 0] = head_to_node_dir.x * lightForce * LightIntensityMutiplier; // To Do: should consider devide by mass
                light_force_array[3 * i + 1] = head_to_node_dir.y * lightForce * LightIntensityMutiplier;
                light_force_array[3 * i + 2] = head_to_node_dir.z * lightForce * LightIntensityMutiplier;
            }
            else {
                light_force_array[3 * i + 0] = 0.0f;
                light_force_array[3 * i + 1] = 0.0f;
                light_force_array[3 * i + 2] = 0.0f;
            }
        }
        ClothSimImp.light_force = Matrix<double>.Build.DenseOfColumnMajor(nHairs * nNodesPerHair * 3, 1, light_force_array);
    }

    void LocalGlobal_update_gravityForce()
    {
        ClothSimImp.gravity_force = gravity * ClothSimImp.gravity_unit;
    }

    void LocalGlobal_update_externalForce()
    {
        LocalGlobal_update_lightForce();
        LocalGlobal_update_gravityForce();
        ClothSimImp.external_force = ClothSimImp.mass_matrix * (ClothSimImp.light_force + ClothSimImp.gravity_force); 
    }

    void LocalGlobal_update_matrixL()
    {
        double sk = Convert.ToDouble(stiffness);
        ClothSimImp.L_matrix = sk * ClothSimImp.ATA_matrix;
    }

    void LocalGlobal_update_matrixJ()
    {
        double sk = Convert.ToDouble(stiffness);
        ClothSimImp.J_matrix = sk * ClothSimImp.AT_matrix;
    }

    Matrix<double> LocalGlobal_evaluate_d()
    {
        Matrix<double> d = Matrix<double>.Build.Dense(nEdges * 3, 1, 0.0f);
        for (int i = 0; i < nEdges; i++)
        {
            SpringConstraint s = ClothSimImp.springConstraint[i];
            Matrix<double> p1 = ClothSimImp.curr_p.SubMatrix(3*s.m_v1, 3, 0, 1); 
            Matrix<double> p2 = ClothSimImp.curr_p.SubMatrix(3*s.m_v2, 3, 0, 1);
            Matrix<double> di = (p1 - p2).NormalizeColumns(2.0) * nodeDistance;
            d.SetSubMatrix(3*i, 3, 0, 1, di);
        }
        return d;
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

    public void updateNodeDistance(float value)
    {
        nodeDistance = value;
        if (ClothSimImp.pinConstraint == null) { return; }
        for (int i = 0; i < ClothSimImp.pinConstraint.Length; i++) // update the pinned node position.
        {
            int nodeIndex = ClothSimImp.pinConstraint[i].m_v0;
            int row = nodeIndex % nHairs; // which row
            int col = nodeIndex / nHairs; // whichcolumn

            double[] new_pos = {  nodeDistance * (col - nHairs / 2),
                               -nodeDistance * (row - nNodesPerHair / 2),
                                0.0f};
            Matrix<double> new_pos_matrix = Matrix<double>.Build.DenseOfColumnMajor(3, 1, new_pos);

            // replace position sub matrix:
            ClothSimImp.rest_p.SetSubMatrix(3 * nodeIndex, 3, 0, 1, new_pos_matrix);
        }
    }

    public void updateLightLookAtSettings(Vector3 lightHeadPos, Vector3 lightHeadLootAt)
    {
        headPos = lightHeadPos;
        headLookAtPos = lightHeadLootAt;
    }

    public void updatelightForce(float f)
    {
        lightForce = f;
    }

    public void updatelightAngle(float a)
    {
        lightAngle = a;
    }

    Vector3 getVector3FromColumnMatrix(Matrix<double> columnMatrix, int index)
    {
        //columnMatrix must be a 3m x 1 matrix.
        //Debug.Assert(columnMatrix.ColumnCount == 1);
        //Debug.Assert(columnMatrix.RowCount == 3 * clothVertArray.Length);
        return new Vector3(Convert.ToSingle(columnMatrix.At(3 * index, 0)),
                           Convert.ToSingle(columnMatrix.At(3 * index + 1, 0)),
                           Convert.ToSingle(columnMatrix.At(3 * index + 2, 0)));
    }

    void setColumnMatrixFromVector3(Matrix<double> columnMatrix, int index, Vector3 vec3)
    {
        //columnMatrix must be a 3m x 1 matrix.
        //Debug.Assert(columnMatrix.ColumnCount == 1);
        //Debug.Assert(columnMatrix.RowCount == 3 * clothVertArray.Length); 
        Matrix<double> new_sub_p_matrix = Matrix<double>.Build.Dense(3, 1, new double[] { vec3.x, vec3.y, vec3.z});
        columnMatrix.SetSubMatrix(3 * index, 3, 0, 1, new_sub_p_matrix);
        columnMatrix.SetSubMatrix(3 * index, 3, 0, 1, new_sub_p_matrix);
        columnMatrix.SetSubMatrix(3 * index, 3, 0, 1, new_sub_p_matrix);
    }

    //////////////////////////// Not used functions ////////////////////////////
    /*
     * examples of using the Math library.
     * Use this to test the math net functions.
     */
    void codeMatrixExample()
    {
        // void SetSubMatrix(int rowIndex, int rowCount, int columnIndex, int columnCount, Matrix<T> subMatrix)
        // Matrix<double> p1 = ClothSimImp.curr_p.SubMatrix(3*s.m_v1, 3, 0, 1); 

        // void SetSubMatrix(int rowIndex, int rowCount, int columnIndex, int columnCount, Matrix<T> subMatrix)
        // d.SetSubMatrix(3*i, 3, 0, 1, di);

        // basic operation
        Debug.Log("Matrix Addition:");
        Matrix<double> A = DenseMatrix.OfArray(new double[,] {  {1,1},
                                                                {1,2} });
        Matrix<double> B = DenseMatrix.OfArray(new double[,] {  {-1,-1},
                                                                {3, 4} });
        Debug.Log("A + B = " + (A + B));

        // build Matrix from Vector
        Debug.Log(" Build Matrix from Vector { 1, 2, 3, 4 }: ");
        double[] v = { 1, 2, 3, 4 };
        Matrix<double> C = Matrix<double>.Build.DenseOfColumnMajor(2, 2, v);
        Matrix<double> D = Matrix<double>.Build.DenseOfColumnMajor(4, 1, v);
        Matrix<double> E = Matrix<double>.Build.DenseOfColumnMajor(1, 4, v);
        Debug.Log(C);
        Debug.Log(D);
        Debug.Log(E);

        // build from index
        Debug.Log("Build Matrix from row/col index: (0,0,2), (0,1,-3)");
        Tuple<int, int, double>[] t = { Tuple.Create(0, 0, 2.0), Tuple.Create(0, 1, -3.0) };
        Matrix<double> F = Matrix<double>.Build.DenseOfIndexed(3, 4, t);
        Debug.Log(F);

        // Matrix Operation
        Debug.Log("Matrix Operation");
        Matrix<double> DE = D * E;
        Matrix<double> FD = F * D;
        Debug.Log("D * E = " + DE);
        Debug.Log("F * D = " + FD);

        // Modify the elements of a Matrix
        Debug.Log("Replace Matrix elements");
        Matrix<double> SUB = Matrix<double>.Build.DenseOfArray(new double[,] { { 1983 } }); // 1 element matrix
        Debug.Log("SUB :" + SUB);
        C.SetSubMatrix(0, 1, 0, 1, SUB); // void SetSubMatrix(int rowIndex, int rowCount, int columnIndex, int columnCount, Matrix<T> subMatrix)
        Debug.Log("C:" + C);

        // solve equation Hx=b
        var H = Matrix<double>.Build.DenseOfArray(new double[,] {   { 3, 2, -1 },
                                                                    { 2, -2, 4 },
                                                                    { -1, 0.5, -1 } });
        var b = Vector<double>.Build.Dense(new double[] { 1, -2, 0 });
        var b2 = Matrix<double>.Build.DenseOfArray(new double[,] {  { 1 },
                                                                    { -2 },
                                                                    { 0 } });
        var x = H.Solve(b2);
        Debug.Log("x:" + x);
    }
}
