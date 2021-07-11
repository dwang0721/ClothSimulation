using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SmallCameraBehavior : MonoBehaviour
{
    public LineBehavior cameraForward;
    public SceneNode attachSceneNode;

    // Start is called before the first frame update
    void Start()
    {
        Debug.Assert(cameraForward);
        Debug.Assert(attachSceneNode);
    }

    // Update is called once per frame
    void Update()
    {
        lookforward();
    }

    void lookforward()
    {
        // draw line
        Matrix4x4 xform = attachSceneNode.getCombinedParentXform();
        Vector3 start = xform.MultiplyPoint(cameraForward.initialPosition);
        Vector3 end = -5.0f * attachSceneNode.transform.forward + start;
        cameraForward.lineUpdate(start, end, new Color(0.7f, 0.0f, 0.5f, 1.0f), 0.1f);

        // place camera
        gameObject.transform.position = start;
        MyLookRotation(Vector3.Normalize(end-start), attachSceneNode.transform.up);
    }

    void MyLookRotation(Vector3 w, Vector3 v)
    {
        float angle = Mathf.Acos(Vector3.Dot(Vector3.up, v)) * Mathf.Rad2Deg;
        Vector3 axis = Vector3.Cross(Vector3.up, v);
        transform.rotation = Quaternion.AngleAxis(angle, axis);

        angle = Mathf.Acos(Vector3.Dot(gameObject.transform.forward, w)) * Mathf.Rad2Deg;
        axis = Vector3.Cross(gameObject.transform.forward, w);
        transform.rotation = Quaternion.AngleAxis(angle, axis) * transform.rotation;
    }
}
