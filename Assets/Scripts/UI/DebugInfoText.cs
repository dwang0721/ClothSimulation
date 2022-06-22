// This script defines the behavior of the debug information.
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.IO;

public class DebugInfoText : MonoBehaviour
{
    Text debugText;
    GlobalData globalData;

    string staticInfo = "";
    string runtimeInfo = "";
    double timeSum = 0;
    int frameCount = 0;

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
        updateDebugInfoOnUI();

        //Debug.Log("globalData.prifileMode " + globalData.profileMode);
        if (globalData.profileMode) 
        {
            logDebugInfoToFile();
        }
    }

    void updateDebugInfoOnUI()
    {
        double deltaTime = Time.unscaledDeltaTime;
        double fps = 1.0f / Time.unscaledDeltaTime;

        runtimeInfo = "";
        runtimeInfo += "[Delta Time] \t" + deltaTime.ToString("F2") + "\n";
        runtimeInfo += "[FPS] \t" + fps.ToString("F2") + "\n";

        debugText.text = staticInfo + runtimeInfo;
    }

    void logDebugInfoToFile() 
    {
        if (frameCount < globalData.profileFrameMax)
        {
            frameCount++;
            if (frameCount >= globalData.profileFrameMin)
            {
                timeSum += Time.unscaledDeltaTime;
            }
        }

        if (frameCount == globalData.profileFrameMax)
        {
            double deltaTimeAvg = timeSum / globalData.profileFrameMax;
            double fpsAvg = 1.0f / deltaTimeAvg;

            string logContent = "";
            logContent += globalData.simulationMode + "\t";
            logContent += globalData.resolution + " x " + globalData.resolution + " : \t";
            logContent += "Delta: " + deltaTimeAvg + ", \t" + "FPS: " + fpsAvg + "\n";
            Debug.Log("Write to file:" + logContent + " at frame " + frameCount);
            WriteDebugInfoToFile(logContent);
            frameCount++;
        }
    }

    void WriteDebugInfoToFile(string content)
    {
        // path of the file
        string path = Application.dataPath + "/Log.txt";

        //  create file if does not exist
        if (!File.Exists(path))
        {
            File.WriteAllText(path, "Peformance Metrics \nLog date: " + System.DateTime.Now + "\n\n");
        }

        File.AppendAllText(path, content);
    }
}
