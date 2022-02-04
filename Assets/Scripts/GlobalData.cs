using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public enum SimulationMethod { ExplicitGPU, LocalGlobal };

/*
 * This is the static simulation controller. The simulation properties should not be updated each frame.
 * Data needs to be preserved acrossed different simulation scene:
 * 1. cloth resolution.
 * 2. Simulation method.
 */

public class GlobalData : MonoBehaviour
{
    // Start is called before the first frame update
    public int resolution;
    public SimulationMethod simulationMode; // default to GPU
   
    void Awake()
    {
        resolution = 10;
        simulationMode = SimulationMethod.LocalGlobal;
        Debug.Log("sim:" + simulationMode.ToString());

        GameObject[] objs = GameObject.FindGameObjectsWithTag("GlobalData");

        if (objs.Length > 1) {
            Destroy(this.gameObject);
        }

        DontDestroyOnLoad(gameObject);
    }

}
