using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

public class CurveInterpolator : MonoBehaviour
{
    public Vector3[] keyPoints;
    public Color keyPointColor, interpolateColor;
    static GameObject[] keyPointGeos;
    public GameObject keyPointPrefab;
    static GameObject interpolatedPoint;

    [Range(0, 0.999f)]
    public float t;

    public bool randomGeneration = false;
    public float nodeDistance = 3.0f;

    // catmull-Rom Matrix
    static Matrix<double> catmullRomMatrix = 0.5 * Matrix<double>.Build.DenseOfArray(
        new double[,] {     { -1,  3, -3,  1 },
                            { 2,  -5,  4, -1 },
                            { -1,  0,  1,  0 },
                            {  0,  2,  0,  0 }, });
    
    // 
    private Matrix<double> U = Matrix<double>.Build.DenseOfArray(new double[,] { { 0, 0, 0, 0 } });

    // Start is called before the first frame update
    void Start()
    {
        Debug.Assert(keyPointPrefab);

        populateKeyPoints(10);
        Debug.Assert(keyPoints.Length >= 4);
        
        initKeyPointGeos();
    }

    // Update is called once per frame
    void Update()
    {
        updateKeyPointGeo();
        catmullRomInterpolation(t);
    }

    void populateKeyPoints(int numKeyPoints)
    {
        
        if (keyPoints.Length < 4 || randomGeneration) {
            keyPoints = new Vector3[numKeyPoints];

            for (int i = 0; i < numKeyPoints; i++) {
                float x = Random.Range(-0.1f, 0.1f) + nodeDistance * (i - numKeyPoints / 2.0f);
                float y = Random.Range(-5.0f, 5.0f);
                float z = -30.0f + Random.Range(-2.0f, 2.0f);
                keyPoints[i] = new Vector3(x, y, z);
            }

        }
    }

    void initKeyPointGeos()
    {
        keyPointGeos = new GameObject[keyPoints.Length];
        for (int i = 0; i < keyPoints.Length; i++) {
            Vector3 location = keyPoints[i];
            var newGeo = Instantiate(keyPointPrefab, location, Quaternion.identity);
            newGeo.transform.localScale = new Vector3(1, 1, 1) * 0.8f;
            newGeo.GetComponent<Renderer>().material.color = keyPointColor;
            keyPointGeos[i] = newGeo;
        }

        interpolatedPoint = Instantiate(keyPointPrefab, keyPoints[0], Quaternion.identity);
        interpolatedPoint.GetComponent<Renderer>().material.color = interpolateColor;
        interpolatedPoint.transform.localScale = new Vector3(1, 1, 1) * 1.0f;
    }

    void updateKeyPointGeo()
    {
        for (int i = 0; i < keyPoints.Length; i++) {
            keyPointGeos[i].transform.position = keyPoints[i];
        }
    }

    void catmullRomInterpolation(float t)
    {
        // t is the global interpolation
        int numSegment = keyPoints.Length - 3; // the number of curve segments, do not count the front and end points.
        float u = (t % (1.0f / numSegment)) / (1.0f / numSegment); // the interpolation value in segmented curve.
        int curveIndex = (int) Mathf.Floor(t / (1.0f / numSegment) + 1);

        U = Matrix<double>.Build.DenseOfArray(new double[,] { { u * u * u, u * u, u, 1 } });

        Matrix<double> keyPointX = Matrix<double>.Build.DenseOfArray(new double[,] {    { keyPoints[curveIndex - 1].x },
                                                                                        { keyPoints[curveIndex].x },
                                                                                        { keyPoints[curveIndex + 1].x }, 
                                                                                        { keyPoints[curveIndex + 2].x } });

        Matrix<double> keyPointY = Matrix<double>.Build.DenseOfArray(new double[,] {    { keyPoints[curveIndex - 1].y },
                                                                                        { keyPoints[curveIndex].y },
                                                                                        { keyPoints[curveIndex + 1].y },
                                                                                        { keyPoints[curveIndex + 2].y } });

        Matrix<double> keyPointZ = Matrix<double>.Build.DenseOfArray(new double[,] {    { keyPoints[curveIndex - 1].z },
                                                                                        { keyPoints[curveIndex].z },
                                                                                        { keyPoints[curveIndex + 1].z },
                                                                                        { keyPoints[curveIndex + 2].z } });


        float x = (float)(U * catmullRomMatrix * keyPointX).At(0, 0);
        float y = (float)(U * catmullRomMatrix * keyPointY).At(0, 0);
        float z = (float)(U * catmullRomMatrix * keyPointZ).At(0, 0);

        interpolatedPoint.gameObject.transform.position = new Vector3(x, y, z);
    }
}
