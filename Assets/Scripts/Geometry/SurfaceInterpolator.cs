using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class SurfaceInterpolator : MonoBehaviour
{
    public Vector3[] keyPoints, verticalMidPoints, horizontalMidPoints;
    public float distanceToCloth = -30.0f;
    public Color keyPointColor, interpolateColor;
    static GameObject[] surfaceKeyPointGeos, verticalMidPointsGeos, horizontalMidPointGeos;
    public GameObject keyPointPrefab;

    bool ready = false;

    [Range(0, 0.999f)]
    public float tx, ty;

    int numX = CPU3D.nHairs - 3;
    int numY = CPU3D.nNodesPerHair - 3;

    private Matrix<double> catmullRomMatrix = 0.5 * Matrix<double>.Build.DenseOfArray(
        new double[,] {     { -1,  3, -3,  1 },
                            { 2,  -5,  4, -1 },
                            { -1,  0,  1,  0 },
                            {  0,  2,  0,  0 }, });

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        populateSurfaceKeyPoints();
    }

    void populateSurfaceKeyPoints()
    {
        if (CPU3D.hairNodesArray.Length == 0) {
            return;
        }


        if (ready) {
            return;
        }

        surfaceKeyPointGeos = new GameObject[CPU3D.nHairs * CPU3D.nNodesPerHair];
        keyPoints = new Vector3[CPU3D.nHairs * CPU3D.nNodesPerHair];

        for (int i = 0; i < surfaceKeyPointGeos.Length; i++)
        {
            Vector3 location = new Vector3(CPU3D.hairNodesArray[i].x, CPU3D.hairNodesArray[i].y, CPU3D.hairNodesArray[i].z + distanceToCloth);
            keyPoints[i] = location;
            var newitem = Instantiate(keyPointPrefab, location, Quaternion.identity);
            newitem.transform.localScale = new Vector3(1, 1, 1) * 0.3f;
            newitem.GetComponent<Renderer>().material.color = keyPointColor;
            surfaceKeyPointGeos[i] = newitem;
            
        }

        // populateInterpolatePoints();

        ready = true;
    }

    //void populateInterpolatePoints()
    //{
    //    interplolatePointsGeos = new GameObject[numX * numY];

    //    for (int i = 0; i < interplolatePointsGeos.Length; i++)
    //    {
    //        Vector3 location = new Vector3(keyPoints[0].x, keyPoints[0].y, keyPoints[0].z);
    //        var newitem = Instantiate(keyPointPrefab, location, Quaternion.identity);
    //        newitem.transform.localScale = new Vector3(1, 1, 1) * 0.5f;
    //        newitem.GetComponent<Renderer>().material.color = interpolateColor;
    //        interplolatePointsGeos[i] = newitem;
    //    }
    //}

    void updateInterpolationPoints()
    {
        for (int i = 0; i < numX; i++) {
            for (int j = 0; j < numY; j++) {
                float tx = (i * 1.0f + 0.5f) / numX;
                float ty = (j * 1.0f + 0.5f) / numY;
                catmullRomInterpolation2D(tx, ty);
            }
        }
    }

    void catmullRomInterpolation2D(float tx, float ty)
    {
        // t is the global interpolation
        int numSegmentX = CPU3D.nHairs - 3;
        float lengthSegmenX = 1.0f / numSegmentX;
        float ux = (tx % lengthSegmenX) / lengthSegmenX; 
        int curveIndexInX = (int)Mathf.Floor(tx / lengthSegmenX + 1);


        int numSegmentY = CPU3D.nHairs - 3;
        float lengthSegmenY = 1.0f / numSegmentY;
        float uy = (ty % lengthSegmenY) / lengthSegmenY;
        int curveIndexInY = (int)Mathf.Floor(ty / lengthSegmenY + 1);

    }


    Vector3 evaluate(float u, Vector3 P0, Vector3 P1, Vector3 P2, Vector3 P3 )
    {
        Matrix<double> U = Matrix<double>.Build.DenseOfArray(new double[,] { { u * u * u, u * u, u, 1 } });

        Matrix<double> keyPointX = Matrix<double>.Build.DenseOfArray(new double[,] {    { P0.x },
                                                                                        { P1.x },
                                                                                        { P2.x },
                                                                                        { P3.x } });

        Matrix<double> keyPointY = Matrix<double>.Build.DenseOfArray(new double[,] {    { P0.y },
                                                                                        { P1.y },
                                                                                        { P2.y },
                                                                                        { P3.y } });

        Matrix<double> keyPointZ = Matrix<double>.Build.DenseOfArray(new double[,] {    { P0.z },
                                                                                        { P1.z },
                                                                                        { P2.z },
                                                                                        { P3.z } });


        float x = (float)(U * catmullRomMatrix * keyPointX).At(0, 0);
        float y = (float)(U * catmullRomMatrix * keyPointY).At(0, 0);
        float z = (float)(U * catmullRomMatrix * keyPointZ).At(0, 0);

        return new Vector3(x, y, z);
    }
}
