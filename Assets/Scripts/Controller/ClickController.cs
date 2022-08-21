using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ClickController : MonoBehaviour
{

    // Update is called once per frame
    void Update()
    {
        // Not in camera rotation mode, and the collider is clicked.
        if (!Input.GetKey(KeyCode.LeftAlt) && Input.GetMouseButtonDown(0))
        {
            RaycastHit hit;
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

            if (Physics.Raycast(ray, out hit))
            {
                GameObject hitObj = hit.collider.gameObject;
                if (hitObj.tag == "Collider")
                {
                    colliderController collider = hitObj.GetComponent<colliderController>();
                    if (collider != null) {
                        collider.setColliderStatusOnClick();
                    }
                }

                if (hitObj.tag == "VertexNode")
                {
                    NodeController node = hitObj.GetComponent<NodeController>();
                    if (node != null) {
                        node.setPinOnClick();
                    }

                }
            }
        }

    }
}
