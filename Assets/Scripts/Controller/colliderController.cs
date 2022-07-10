using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class colliderController : MonoBehaviour
{
    bool isSelected;
    Vector3 offset;
    Vector3 screenPoint;
    static float minSize = 5.0f;
    static float maxSize = 25.0f;

    // Start is called before the first frame update
    void Start()
    {
        gameObject.GetComponent<Renderer>().material.color = new Color(0.5f, 0, 0, 1);
    }

    // Update is called once per frame
    void Update()
    {
        //// Not in camera rotation mode, and the collider is clicked.
        //if (!Input.GetKey(KeyCode.LeftAlt) && Input.GetMouseButtonDown(0)) {
        //    RaycastHit hit;
        //    Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

        //    if (Physics.Raycast(ray, out hit)) {
        //        GameObject hitObj = hit.collider.gameObject;
        //        if (hitObj.tag == "Collider" && hitObj == gameObject) {
        //            screenPoint = Camera.main.WorldToScreenPoint(transform.position);
        //            offset = transform.position - Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z)); ;
        //            isSelected = true;
        //            gameObject.GetComponent<Renderer>().material.color = new Color(0.8f, 0, 0, 1);
        //        }
        //    }
        //}

        if (Input.GetMouseButtonUp(0)) {
            isSelected = false;
            gameObject.GetComponent<Renderer>().material.color = new Color(0.5f, 0, 0, 1);
        }

        updatePosition();
    }

    public void setColliderStatusOnClick()
    {
        screenPoint = Camera.main.WorldToScreenPoint(transform.position);
        offset = transform.position - Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z));
        isSelected = true;
        gameObject.GetComponent<Renderer>().material.color = new Color(0.8f, 0, 0, 1);
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

        if (Input.mouseScrollDelta.y != 0.0f && transform.localScale.x + Input.mouseScrollDelta.y > minSize && transform.localScale.x + Input.mouseScrollDelta.y < maxSize)
        {
            Debug.Assert(transform.localScale.x == transform.localScale.y);
            Debug.Assert(transform.localScale.y == transform.localScale.z);
            transform.localScale = new Vector3( Input.mouseScrollDelta.y + transform.localScale.x, Input.mouseScrollDelta.y + transform.localScale.y, Input.mouseScrollDelta.y + transform.localScale.z);
        }

        Vector3 currentScreenPoint = new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z);
        Vector3 currentPosition = Camera.main.ScreenToWorldPoint(currentScreenPoint) + offset;
        transform.position = currentPosition;
    }
}
