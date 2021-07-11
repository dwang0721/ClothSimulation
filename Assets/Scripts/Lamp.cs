using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Lamp : MonoBehaviour
{
    public CPU3D simulation;
    public float tumbleSpeed;
    float viewingAngle;
    public Material lightMaterial;
    public Transform head;
    public Vector3 headLookAt;
    public GameObject lineSegmentPrefab;

    // Line Segments
    private GameObject lsCenter;
    private GameObject lsLeft;
    private GameObject lsRight;
    private GameObject lsUp;
    private GameObject lsDown;

    bool gotCollider;
    bool isSelected;
    bool headSelected;
    Vector3 offset;
    Vector3 screenPoint;
    GameObject collider;
    

    // Start is called before the first frame update
    void Start()
    {        
        Debug.Assert(simulation);
        Debug.Assert(head);        

        // Instantiate line segments
        lsCenter = Instantiate(lineSegmentPrefab, Vector3.zero, Quaternion.identity);
        lsRight = Instantiate(lineSegmentPrefab, Vector3.zero, Quaternion.identity);
        lsLeft = Instantiate(lineSegmentPrefab, Vector3.zero, Quaternion.identity);
        lsUp = Instantiate(lineSegmentPrefab, Vector3.zero, Quaternion.identity);
        lsDown = Instantiate(lineSegmentPrefab, Vector3.zero, Quaternion.identity);
    }

    void updatePosition() {
        if (!isSelected) { return; }

        if (Input.GetKey(KeyCode.LeftControl)) {
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

    void drawSpotLight()
    {
        float angleCos = Vector3.Dot(head.transform.forward, Vector3.forward);
        float headToScreen = Mathf.Abs(0.0f - head.transform.position.z);
        float headToIntersetion = headToScreen / angleCos;
        float lightRadius = headToScreen * Mathf.Tan(viewingAngle * 0.5f * Mathf.Deg2Rad);

        headLookAt = head.transform.position + headToIntersetion * head.transform.forward;
        Vector3 lgithLeft = headLookAt - head.transform.right * lightRadius;
        Vector3 lgithRight = headLookAt + head.transform.right * lightRadius;
        Vector3 lgithUp = headLookAt + head.transform.up * lightRadius;
        Vector3 lgithDown = headLookAt - head.transform.up * lightRadius;

        updateLineSegment(lsCenter, head.transform.position, headLookAt);
        updateLineSegment(lsLeft, head.transform.position, lgithLeft);
        updateLineSegment(lsRight, head.transform.position, lgithRight);
        updateLineSegment(lsUp, head.transform.position, lgithUp);
        updateLineSegment(lsDown, head.transform.position, lgithDown);
    }

    void updateLineSegment(GameObject ls, Vector3 start, Vector3 end) {
        Vector3 directionVector = end - start;
        Vector3 scale = new Vector3(ls.transform.localScale.x, directionVector.magnitude / 2, ls.transform.localScale.z);
        ls.transform.localScale = scale;
        ls.transform.rotation = Quaternion.FromToRotation(Vector3.up, directionVector.normalized);
        ls.transform.position = start + ls.transform.up * ls.transform.localScale.y;
    }

    void updateSpotLight()
    {
        Light spotLight = GameObject.Find("Spot Light").GetComponent<Light>();
        spotLight.spotAngle = viewingAngle;
    }

    void tumbleHead(float xDirection, float yDirection) {
        // store z component so we can lock the z-axis
        
        float zComp = head.transform.localRotation.eulerAngles.z;
        Quaternion yaw = Quaternion.AngleAxis(xDirection * tumbleSpeed, head.transform.up);
        Quaternion pitch = Quaternion.AngleAxis(yDirection * tumbleSpeed, head.transform.right);

        Quaternion q = pitch * yaw;
        Quaternion newRot = q * head.transform.localRotation;
        Vector3 newEuler = newRot.eulerAngles;
        newEuler.z = zComp;
        newRot.eulerAngles = newEuler;
        head.transform.localRotation = newRot;
    }

    void updateHeadRotation() {
        if (!headSelected) { return; }

        tumbleHead(Input.GetAxis("Mouse X"), -Input.GetAxis("Mouse Y"));
    }

    void updateLightSimulation()
    {
        simulation.updateLightDirection(head.transform.position, headLookAt);
    }

    public void setViewingAngle(float angle)
    {
        viewingAngle = angle;
    }

    public void setLightColor(float intensity)
    {
        Light spotLight = GameObject.Find("Spot Light").GetComponent<Light>();
        float lerp = intensity / 0.3f ;
        spotLight.color = Color.red* lerp + Color.green * (1- lerp);
    }

    // Update is called once per frame
    void Update()
    {   
        if (!gotCollider) {
            collider = simulation.getCollider();
            if (collider != null) {
                gotCollider = true;
            }
        }

        if (!Input.GetKey(KeyCode.LeftAlt) && Input.GetMouseButtonDown(0)) {
            RaycastHit hit;
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

            if (Physics.Raycast(ray, out hit)) {
                GameObject hitObj = hit.collider.gameObject;
                if (hitObj.name == "base" && hitObj.GetComponentInParent<Lamp>() != null) {
                    screenPoint = Camera.main.WorldToScreenPoint(transform.position);
                    offset = transform.position - Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z)); ;
                    isSelected = true;

                }

                if (hitObj.name == "head") {
                    headSelected = true;
                }
            }
        }

        if (Input.GetMouseButtonUp(0)) {
            isSelected = false;
            headSelected = false;
        }

        updatePosition();
        updateHeadRotation();
        updateSpotLight();
        updateLightSimulation();
        drawSpotLight();        
    }
}
