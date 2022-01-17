using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class colliderController : MonoBehaviour
{
    bool isSelected;
    Vector3 offset;
    Vector3 screenPoint;

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        if (!Input.GetKey(KeyCode.LeftAlt) && Input.GetMouseButtonDown(0)) {
            RaycastHit hit;
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

            if (Physics.Raycast(ray, out hit)) {
                GameObject hitObj = hit.collider.gameObject;
                if (hitObj.tag == "Collider") {
                    screenPoint = Camera.main.WorldToScreenPoint(transform.position);
                    offset = transform.position - Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z)); ;
                    isSelected = true;
                }
            }
        }

        if (Input.GetMouseButtonUp(0)) {
            isSelected = false;
        }

        updatePosition();
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
}
