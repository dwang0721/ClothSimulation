using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FlagController : MonoBehaviour
{
    public CPU3D explictSimulator;
    public SimulationMethod simulationMode;
    public float tumbleSpeed;
    private Transform flagPole;
    private Transform flagCube;

    bool isSelected;
    bool flagPoleSelected;

    Vector3 offset;
    Vector3 screenPoint;

    // Start is called before the first frame update
    void Start()
    {
        flagCube = gameObject.transform;
        flagPole = flagCube.Find("FlagPoleMesh");
        Debug.Assert(explictSimulator);
        Debug.Assert(flagPole);
        Debug.Assert(flagCube);
    }

    // Update is called once per frame
    void Update()
    {
        setClickedStatus();
        updatePosition();
        updateFlagPoleRotation();
    }

    void setClickedStatus()
    {
        if (!Input.GetKey(KeyCode.LeftAlt) && Input.GetMouseButtonDown(0))
        {
            RaycastHit hit;
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

            if (Physics.Raycast(ray, out hit))
            {
                GameObject hitObj = hit.collider.gameObject;
                if (hitObj.name == "FlagControl" && hitObj.GetComponentInParent<FlagController>() != null)
                {
                    screenPoint = Camera.main.WorldToScreenPoint(transform.position);
                    offset = transform.position - Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z)); ;
                    isSelected = true;
                }

                if (hitObj.name == "FlagPoleMesh")
                {
                    flagPoleSelected = true;
                }
            }
        }

        if (Input.GetMouseButtonUp(0))
        {
            isSelected = false;
            flagPoleSelected = false;
        }
    }

    void updatePosition()
    {
        if (!isSelected) { return; }

        if (Input.GetKey(KeyCode.LeftControl))
        {
            // Handle z-axis
            float xDirection = Input.GetAxis("Mouse X");
            Vector3 mv = new Vector3(0f, 0f, xDirection * 100 * Time.deltaTime);
            transform.position += mv;
            return;
        }

        Vector3 currentScreenPoint = new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z);
        Vector3 currentPosition = Camera.main.ScreenToWorldPoint(currentScreenPoint) + offset;
        transform.position = currentPosition;
    }

    void tumbleFlagPole(float xDirection, float yDirection)
    {
        // store x component so we can lock the x-axis
        float xComp = flagPole.transform.localRotation.eulerAngles.x;
        Quaternion yaw = Quaternion.AngleAxis(xDirection * tumbleSpeed, flagPole.transform.up); //y
        Quaternion roll = Quaternion.AngleAxis(yDirection * tumbleSpeed, flagPole.transform.forward); // z

        Quaternion q = roll * yaw;
        Quaternion newRot = q * flagPole.transform.localRotation;
        Vector3 newEuler = newRot.eulerAngles;
        newEuler.x = xComp;
        newRot.eulerAngles = newEuler;
        flagPole.transform.localRotation = newRot;
    }

    void updateFlagPoleRotation()
    {
        if (!flagPoleSelected) { return; }
        tumbleFlagPole(-Input.GetAxis("Mouse X"), Input.GetAxis("Mouse Y"));
    }

    public void setGlobalPose(Vector3 pos)
    {
        flagCube.position = pos;
    }

    public Quaternion getPoleRotation()
    {
        return flagPole.transform.rotation;
    }
}
