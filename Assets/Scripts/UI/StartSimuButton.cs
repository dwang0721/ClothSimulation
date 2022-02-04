// This script defines the "Start Simulation" button behaviors

using System.Collections;
using System.Collections.Generic;
using UnityEngine.SceneManagement;
using UnityEngine;
using UnityEngine.UI;

public class StartSimuButton : MonoBehaviour
{
    // Start is called before the first frame update
    Button btn;
    GlobalData globalData;
    public Dropdown resolutionDropDown;
    public Dropdown simulatorDropDown;

    void Awake()
    {
        btn = gameObject.GetComponent<Button>();
        btn.onClick.AddListener(delegate { TaskOnClick(); });
    }

    void Start() {
        globalData = GameObject.FindGameObjectsWithTag("GlobalData")[0].GetComponent<GlobalData>();
        Debug.Assert(globalData);

        // load the UI data and update the dropdown selection
        updateResolutionDropDownOnUI();
        updateSimulatornDropDownOnUI();
    }

    // Update is called once per frame
    void TaskOnClick()
    {
        // set the global values
        updateResolutionOnGobal();
        updateSimulatorOnGobal();

        // reload the scene
        SceneManager.LoadScene("RT_Sample");
    }

    void updateResolutionDropDownOnUI() {
        switch (globalData.resolution) 
        {
            case 8:
                resolutionDropDown.value = 0;
                break;
            case 10:
                resolutionDropDown.value = 1;
                break;
            case 12:
                resolutionDropDown.value = 2;
                break;
            case 16:
                resolutionDropDown.value = 3;
                break;
            case 20:
                resolutionDropDown.value = 4;
                break;
            case 24:
                resolutionDropDown.value = 5;
                break;
            case 32:
                resolutionDropDown.value = 6;
                break;
            case 64:
                resolutionDropDown.value = 7;
                break;
            default:
                resolutionDropDown.value = 0;
                break;
        }
    }

    void updateSimulatornDropDownOnUI()
    {
        switch (globalData.simulationMode)
        {
            case SimulationMethod.ExplicitGPU:
                simulatorDropDown.value = 0;
                break;
            case SimulationMethod.LocalGlobal:
                simulatorDropDown.value = 1;
                break;
            default:
                simulatorDropDown.value = 0;
                break;
        }
    }

    void updateResolutionOnGobal()
    {
        //Debug.Log("Res: " + resolutionDropDown.value);
        switch (resolutionDropDown.value)
        {
            case 0:
                globalData.resolution = 8;
                break;
            case 1:
                globalData.resolution = 10;
                break;
            case 2:
                globalData.resolution = 12;
                break;
            case 3:
                globalData.resolution = 16;
                break;
            case 4:
                globalData.resolution = 20;
                break;
            case 5:
                globalData.resolution = 24;
                break;
            case 6:
                globalData.resolution = 32;
                break;
            case 7:
                globalData.resolution = 64;
                break;
            default:
                globalData.resolution = 10;
                break;
        }
    }

    void updateSimulatorOnGobal() 
    {
        //Debug.Log("Simulator: " + simulatorDropDown.value);
        switch (simulatorDropDown.value)
        {
            case 0:
                globalData.simulationMode = SimulationMethod.ExplicitGPU;
                break;
            case 1:
                globalData.simulationMode = SimulationMethod.LocalGlobal;
                break;
            default:
                globalData.simulationMode = SimulationMethod.LocalGlobal;
                break;
        }
    }
}
