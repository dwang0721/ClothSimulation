using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraBehavior : MonoBehaviour
{
    public Transform targetTransform = null;
    private Vector3 defaultTargetPos, defaultCameraPos;
    private Quaternion defaultCameraRot;
    // Start is called before the first frame update
    void Start()
    {
        Debug.Assert(targetTransform);
        defaultTargetPos = targetTransform.position;
        defaultCameraPos = gameObject.transform.position;
        defaultCameraRot = gameObject.transform.localRotation;
    }

    // Update is called once per frame
    void Update()
    {
        lookAtTarget();
        cameraTumble();
        cameraTrack();
        cameraDolly();
        resetCamera();
    }

    void cameraTumble()
    {
        if (Input.GetMouseButton(0) && Input.GetKey(KeyCode.LeftAlt))
        {
            float x = Input.GetAxis("Mouse X");
            float y = Input.GetAxis("Mouse Y");


            // prevent camera flip over
            Vector3 cameraToTarget = transform.position - targetTransform.position;
            float cameraDotWorldUp = Vector3.Dot(Vector3.up, Vector3.Normalize(cameraToTarget));
            float threshold = 0.85f;
            if (cameraDotWorldUp >= threshold && y < 0.0f) return;
            if (cameraDotWorldUp <= -threshold && y > 0.0f) return;
            Tumble(-y, transform.right);

            Tumble(x, transform.up);    
            
        }
    }

    void cameraTrack()
    {
        if (Input.GetMouseButton(1) && Input.GetKey(KeyCode.LeftAlt))
        {
            float x = Input.GetAxis("Mouse X");
            float y = Input.GetAxis("Mouse Y");

            Vector3 shift = new Vector3(x, y, 0.0f);
            gameObject.transform.position += shift;
            targetTransform.position += shift;
        }
    }

    void cameraDolly()
    {
        if (Input.mouseScrollDelta.y != 0.0f && Input.GetKey(KeyCode.LeftAlt))
        {
            Vector3 targetToCamera = Vector3.Normalize(targetTransform.position - transform.position);
            gameObject.transform.position += targetToCamera * Input.mouseScrollDelta.y;
        }
    }


    void resetCamera()
    {
        if (Input.GetKey(KeyCode.R))
        {
            targetTransform.position = defaultTargetPos;
            gameObject.transform.position = new Vector3 (defaultCameraPos.x, defaultCameraPos.y, defaultCameraPos.z);
            gameObject.transform.localRotation = defaultCameraRot;
            lookAtTarget();
        }
    }

    void Tumble(float degree, Vector3 axis)
    {
        Vector3 cameraToTarget = transform.position - targetTransform.position;

        Quaternion q = Quaternion.AngleAxis(degree, axis);

        Matrix4x4 rot = Matrix4x4.TRS(Vector3.zero, q, Vector3.one);
        Matrix4x4 invP = Matrix4x4.TRS(cameraToTarget, Quaternion.identity, Vector3.one);
        rot = invP.inverse * rot * invP;
        transform.localPosition = rot.MultiplyPoint(transform.localPosition);
        lookAtTarget();
    }

    void lookAtTarget()
    {
        Vector3 w = Vector3.Normalize(targetTransform.position - transform.position);
        Vector3 u = Vector3.Normalize(Vector3.Cross(-w, Vector3.up));
        Vector3 v = Vector3.Normalize(Vector3.Cross(u, -w));
        MyLookRotation(w, v);
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
