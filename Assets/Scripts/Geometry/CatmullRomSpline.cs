using System.Collections;
using System.Collections.Generic;
using UnityEngine;


//https://github.com/JPBotelho/Catmull-Rom-Splines
public class CatmullRomSpline : MonoBehaviour
{
    //Struct to keep position, normal and tangent of a spline point
    [System.Serializable]
    public struct CatmullRomPoint
    {
        public Vector3 position;
        public Vector3 tangent;
        public Vector3 normal;

        public CatmullRomPoint(Vector3 position, Vector3 tangent, Vector3 normal)
        {
            this.position = position;
            this.tangent = tangent;
            this.normal = normal;
        }
    }
}
