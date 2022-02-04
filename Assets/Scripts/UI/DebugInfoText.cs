// This script defines the behavior of the debug information.

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class DebugInfoText : MonoBehaviour
{
    Text debugText;
    GlobalData globalData;

    string staticInfo = "";
    string runtimeInfo = "";

    // Start is called before the first frame update
    void Awake()
    {
        debugText = gameObject.GetComponent<Text>();
    }

    void Start() 
    {
        globalData = GameObject.FindGameObjectsWithTag("GlobalData")[0].GetComponent<GlobalData>();
        Debug.Assert(globalData);

        string resolutionText = "" + globalData.resolution;
        string simulatorText = globalData.simulationMode.ToString();

        debugText.text = "";
        staticInfo += SystemInfo.processorType + "\n";
        staticInfo += "CPU Frequency: " + SystemInfo.processorFrequency + " MHz, \t";
        staticInfo += "Thread Count: " + SystemInfo.processorCount + "\n ----------------------------------- \n";
        staticInfo += "[Simulator] \t" + simulatorText + "\n";
        staticInfo += "[Resolution] \t" + resolutionText + " x " + resolutionText + "\n";
    }

    // Update is called once per frame
    void Update()
    {
        string deltaTime = Time.unscaledDeltaTime.ToString();
        string fps = "" + 1.0f / Time.unscaledDeltaTime;

        runtimeInfo = "";
        runtimeInfo += "[Delta Time] \t" + deltaTime + "\n";
        runtimeInfo += "[FPS] \t" + fps + "\n";

        debugText.text = staticInfo + runtimeInfo;
    }
}
