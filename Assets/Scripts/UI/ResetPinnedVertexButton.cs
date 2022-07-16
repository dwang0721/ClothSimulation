using System.Collections;
using System.Collections.Generic;
using UnityEngine.SceneManagement;
using UnityEngine;
using UnityEngine.UI;

public class ResetPinnedVertexButton : MonoBehaviour
{
    // Start is called before the first frame update
    Button btn;
    public CPU3D simulationController3D;

    void Awake()
    {
        btn = gameObject.GetComponent<Button>();
        btn.onClick.AddListener(delegate { TaskOnClick(); });
    }

    // Update is called once per frame
    void Start()
    {
        Debug.Assert(simulationController3D);
    }

    // Update is called once per frame
    void TaskOnClick()
    {
        simulationController3D.resetPinnedNodeDistance(1);
    }
}
