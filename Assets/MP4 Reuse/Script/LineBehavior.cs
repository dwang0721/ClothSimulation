using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LineBehavior : MonoBehaviour
{
    Quaternion defaultRataion = Quaternion.identity;
    public Vector3 initialRotation = new Vector3(0.0f, 90.0f, -90.0f); // roatate the line to align with Z axis;
    public Vector3 initialPosition = new Vector3(0.0f, 0.0f, 0.0f);

    // Start is called before the first frame update
    void Start()
    {
        defaultRataion.eulerAngles = initialRotation;
    }

    public void lineUpdate(Vector3 start, Vector3 end, Color color, float thickness)
    {
        // calculate vector direction, length and rotation angle
        Vector3 direction = Vector3.Normalize(start - end);
        float length = Vector3.Magnitude(start - end);
        Quaternion rotation = zToVectorAngle(direction);

        // update gameObject
        gameObject.transform.localScale = new Vector3(thickness, length, thickness);
        gameObject.transform.rotation = rotation * defaultRataion;
        gameObject.transform.localPosition = (start + end) * 0.5f;
        gameObject.GetComponent<Renderer>().material.SetColor("_Color", color);
    }

    Quaternion zToVectorAngle(Vector3 vector)
    {
        Vector3 axis = Vector3.Cross(new Vector3(0, 0, 1), Vector3.Normalize(vector));
        float cosTheta = Vector3.Dot(Vector3.Normalize(vector), new Vector3(0, 0, 1));
        float angle = Mathf.Acos(cosTheta);
        return Quaternion.AngleAxis(Mathf.Rad2Deg * angle, axis);
    }
}
