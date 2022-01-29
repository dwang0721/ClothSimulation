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
