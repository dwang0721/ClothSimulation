using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PrimitiveRotate : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        gameObject.transform.localEulerAngles = new Vector3(gameObject.transform.localEulerAngles.x, gameObject.transform.localEulerAngles.y, 45.0f * Mathf.Sin(Time.timeSinceLevelLoad));
    }
}
