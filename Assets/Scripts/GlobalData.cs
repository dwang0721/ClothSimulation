using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public enum SimulationMethod 
{   
    ExplicitGPU, 
    LocalGlobal 
};

public enum RenderMode
{
    Node,
    Mesh
}

/*
 * This is the static simulation controller. The simulation properties should not be updated each frame.
 * Data needs to be preserved acrossed different simulation scene:
 * 1. cloth resolution.
 * 2. Simulation method.
 * 3. Performance profiling Metrics
 */

public class GlobalData : MonoBehaviour
{
    // Start is called before the first frame update
    public int resolution = 10;
    public SimulationMethod simulationMode = SimulationMethod.LocalGlobal; // default to GPU
    public RenderMode renderMode = RenderMode.Mesh;
    public float clothDebugNodeSize = 0.3f;
    public bool isFreeFallMode = true;

    // profiling metrics setting
    public bool profileMode = false; // should start measuring the performance
    public int profileFrameMin = 11;
    public int profileFrameMax = 110;
    public int rezMin = 25;
    public int rezMax = 46;

    private int frameCount = 0;
    private int curRes;
    private int numSimulators;
    private int simIndex = 0;

    void Awake()
    {
        //resolution = 10;
        //simulationMode = SimulationMethod.LocalGlobal;

        GameObject[] objs = GameObject.FindGameObjectsWithTag("GlobalData");

        if (objs.Length > 1) {
            Destroy(this.gameObject);
        }

        DontDestroyOnLoad(gameObject);
    }

    void Start()
    {
        initProfileModeConfig();
    }

    void Update()
    {
        if (profileMode)
        {
            startProfileModeLoop();
        }
    }

    void initProfileModeConfig() 
    {
        frameCount = 0;
        curRes = rezMin;
        simIndex = 1;
        numSimulators = System.Enum.GetValues(typeof(SimulationMethod)).Length;
    }

    /*
     * This method will be turned on when profileMode is set to True.
     * The animation loop will go thru all the simulators and resolutions 
     * DebugInfoText.cs captures the performance metrics and writes to log.
     */
    void startProfileModeLoop()
    {
        if (frameCount < profileFrameMax + 3)
        {
            frameCount++;
            return; // 
        }

        simulationMode = (SimulationMethod)simIndex;
        resolution = curRes;

        frameCount = 0; // reset frame count
        SceneManager.LoadScene("RT_Sample"); // reload the scene

        if (curRes <= rezMax)
        {
            curRes++;
        }

        if (curRes > rezMax && simIndex < numSimulators)
        {
            curRes = rezMin;
            simIndex++;
        }

        if (curRes > rezMax && simIndex >= numSimulators) 
        {
            profileMode = false;
        }

    }
}

