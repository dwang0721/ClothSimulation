using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NodeController : MonoBehaviour
{
    GameObject simulator;
    CPU3D simulationController3D;
    OptMethod simulationOptMethodController;
    bool isPinned = false;
    int nodeIndex;

    // Start is called before the first frame update
    void Start()
    {
        querrySimulator();
    }

    // Update is called once per frame
    void Update()
    {
        updatePinColor();
    }

    public void setNodeIndex(int index)
    {
        nodeIndex = index;
    }

    private void updateNodePinStatus()
    {
        // update global node array pin status
        int status = isPinned ? 1 : 0;
        
        // GameObject[] simulator= GameObject.FindGameObjectsWithTag("ExplicitGPUSimulator");

        if (simulationController3D) {
            simulationController3D.setHairNodesPinStatusAtIndex(nodeIndex, status);
        }else {
            simulationOptMethodController.setHairNodesPinStatusAtIndex(nodeIndex, status);
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

    void updatePinColor()
    {
        if (isPinned)
        {
            gameObject.GetComponent<Renderer>().material.color = new Color(0, 1.0f, 0, 1);
        }
        else {
            gameObject.GetComponent<Renderer>().material.color = new Color(1.0f, 1.0f, 1.0f, 1);
        }
    }

    void querrySimulator()
    {
        GlobalData globalData = GameObject.FindGameObjectsWithTag("GlobalData")[0].GetComponent<GlobalData>();
        Debug.Assert(globalData);

        if (globalData.simulationMode == SimulationMethod.ExplicitGPU)
        {
            simulator = GameObject.FindGameObjectsWithTag("ExplicitGPUSimulator")[0];
            simulationController3D = simulator.GetComponent<CPU3D>();
        }
        else {
            simulator = GameObject.FindGameObjectsWithTag("ImplicitCPUSimulator")[0];
            simulationOptMethodController = simulator.GetComponent<OptMethod>();
        }

        Debug.Assert(simulator);
    }


}
