using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NodeController : MonoBehaviour
{
    CPU3D simulationController3D;
    OptMethod simulationOptMethodController;
    bool isPinned = false;
    int nodeIndex;

    // Start is called before the first frame update
    void Start()
    {
        //simulationController3D = GameObject.FindGameObjectsWithTag("ExplicitCPUSimulator")[0].GetComponent<CPU3D>();
        //Debug.Assert(simulationController3D);
    }

    // Update is called once per frame
    void Update()
    {
        updatePinStatus();
    }

    public void setNodeIndex(int index)
    {
        nodeIndex = index;
    }

    private void updateNodePinStatus()
    {
        // update global node array pin status
        int status = isPinned ? 1 : 0;
        
        GameObject[] simulator= GameObject.FindGameObjectsWithTag("ExplicitGPUSimulator");

        if (simulator.Length != 0) {
            // explict
            simulationController3D = simulator[0].GetComponent<CPU3D>();
            simulationController3D.setHairNodesPinStatusAtIndex(nodeIndex, status);
        }else {
            simulator = GameObject.FindGameObjectsWithTag("ImplicitCPUSimulator");
            simulationOptMethodController = simulator[0].GetComponent<OptMethod>();
        }
    }

    public void setPinOnClick()
    {
        isPinned = !isPinned;
        updateNodePinStatus();
    }

    public void setIsPinned(int pinStatus)
    {
        // we kept pinStatus as int because we store this variable on the GPU side in 4 bytes.
        // each GPU thread read 4 bytes. Doing so to acchieve coalesced memory access.
        bool status = pinStatus == 0 ? false : true;
        isPinned = status;
    }

    void updatePinStatus()
    {
        if (isPinned)
        {
            gameObject.GetComponent<Renderer>().material.color = new Color(0, 1.0f, 0, 1);
        }
        else {
            gameObject.GetComponent<Renderer>().material.color = new Color(1.0f, 1.0f, 1.0f, 1);
        }
    }


}
