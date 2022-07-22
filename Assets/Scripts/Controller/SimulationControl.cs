//This is the run time simulation controller. The simulation properties can be updated in each frame

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class SimulationControl : MonoBehaviour
{
    public CPU3D simulationController3D;
    public OptMethod simulationOptMethodController;
    public Lamp theLamp;
    public SimulationMethod simulationMode;
    Slider SliderStiffness, SliderGravity, SliderDistance, SliderMaxTravel, SliderBendingStiffness, SliderFriction;
    Slider LightIntensity, LightAngle;
    Slider FlagTransX, FlagTransY, FlagTransZ, FlagRotX, FlagRotY, FlagRotZ;

    GlobalData globalData;

    // some values
    float lightMaxIntensityImplict = 0.3f;
    float lightMaxIntensityExplicit = 0.5f;

    private void Awake()
    {
        // take data from the global setting
        globalData = GameObject.FindGameObjectsWithTag("GlobalData")[0].GetComponent<GlobalData>();
        Debug.Assert(globalData);
        simulationMode = globalData.simulationMode;

        Debug.Log("New Simlator: " + globalData.simulationMode + ";\t" + globalData.resolution );

        // set up the run time simulator
        switch (simulationMode)
        {
            case SimulationMethod.ExplicitGPU:
                simulationController3D.gameObject.SetActive(true);
                simulationOptMethodController.gameObject.SetActive(false);
                simulationController3D.resolution = globalData.resolution;
                simulationController3D.simulationScenario = globalData.simulationScenario;
                simulationController3D.clothDebugNodeSize = globalData.clothDebugNodeSize;
                theLamp.gameObject.SetActive(true);
                theLamp.simulationMode = SimulationMethod.ExplicitGPU;
                break;
            case SimulationMethod.LocalGlobal:
                simulationController3D.gameObject.SetActive(false);
                simulationOptMethodController.gameObject.SetActive(true);
                simulationOptMethodController.resolution = globalData.resolution;
                simulationOptMethodController.clothDebugNodeSize = globalData.clothDebugNodeSize;
                theLamp.gameObject.SetActive(true);
                theLamp.simulationMode = SimulationMethod.LocalGlobal;
                break;
            default:
                // code block
                break;
        }
    }

    // Start is called before the first frame update
    void Start()
    {
        // find slider controls
        SliderDistance = gameObject.transform.Find("RunTimeConfigPanel").Find("distance").GetComponent<Slider>();
        SliderStiffness = gameObject.transform.Find("RunTimeConfigPanel").Find("stiffness").GetComponent<Slider>();
        SliderGravity = gameObject.transform.Find("RunTimeConfigPanel").Find("gravity").GetComponent<Slider>();        
        SliderMaxTravel = gameObject.transform.Find("RunTimeConfigPanel").Find("maxTravel").GetComponent<Slider>();
        SliderBendingStiffness = gameObject.transform.Find("RunTimeConfigPanel").Find("bendStiff").GetComponent<Slider>();
        SliderFriction = gameObject.transform.Find("RunTimeConfigPanel").Find("friction").GetComponent<Slider>();

        LightIntensity = gameObject.transform.Find("RunTimeConfigPanel").Find("intensity").GetComponent<Slider>();
        LightAngle = gameObject.transform.Find("RunTimeConfigPanel").Find("angle").GetComponent<Slider>();

        FlagTransX = gameObject.transform.Find("FlagControlPanel").Find("FlagTransX").GetComponent<Slider>();
        FlagTransY = gameObject.transform.Find("FlagControlPanel").Find("FlagTransY").GetComponent<Slider>();
        FlagTransZ = gameObject.transform.Find("FlagControlPanel").Find("FlagTransZ").GetComponent<Slider>();
        FlagRotX = gameObject.transform.Find("FlagControlPanel").Find("FlagRotX").GetComponent<Slider>();
        FlagRotY = gameObject.transform.Find("FlagControlPanel").Find("FlagRotY").GetComponent<Slider>();
        FlagRotZ = gameObject.transform.Find("FlagControlPanel").Find("FlagRotZ").GetComponent<Slider>();

        Debug.Assert(SliderDistance && SliderStiffness && SliderGravity && SliderMaxTravel && SliderBendingStiffness && SliderFriction);
        Debug.Assert(LightIntensity && LightAngle);
        Debug.Assert(FlagTransX && FlagTransY && FlagTransZ && FlagRotX && FlagRotY && FlagRotZ);

        // set up callbacks
        switch (simulationMode)
        {
            case SimulationMethod.ExplicitGPU:
                SliderDistance.onValueChanged.AddListener(delegate { simulationController3D.updateNodeDistance(SliderDistance.value); });
                SliderStiffness.onValueChanged.AddListener(delegate { simulationController3D.updateStiffness(SliderStiffness.value); });
                SliderGravity.onValueChanged.AddListener(delegate { simulationController3D.updateGravity(SliderGravity.value); });
                SliderMaxTravel.onValueChanged.AddListener(delegate { simulationController3D.updateMaxTravel(SliderMaxTravel.value); });
                SliderBendingStiffness.onValueChanged.AddListener(delegate { simulationController3D.updateBendStiff(SliderBendingStiffness.value); });
                SliderFriction.onValueChanged.AddListener(delegate { simulationController3D.updateFriction(SliderFriction.value); });

                LightIntensity.onValueChanged.AddListener(delegate { simulationController3D.updatelightForce(LightIntensity.value); });
                LightIntensity.onValueChanged.AddListener(delegate { theLamp.setLightColor(LightIntensity.value, lightMaxIntensityImplict); });
                LightAngle.onValueChanged.AddListener(delegate { simulationController3D.updatelightAngle(LightAngle.value); });
                LightAngle.onValueChanged.AddListener(delegate { theLamp.setViewingAngle(LightAngle.value); });
                break;
            case SimulationMethod.LocalGlobal:
                SliderDistance.onValueChanged.AddListener(delegate { simulationOptMethodController.updateNodeDistance(SliderDistance.value); });
                SliderStiffness.onValueChanged.AddListener(delegate { simulationOptMethodController.updateStiffness(SliderStiffness.value); });
                SliderGravity.onValueChanged.AddListener(delegate { simulationOptMethodController.updateGravity(SliderGravity.value); });

                LightIntensity.onValueChanged.AddListener(delegate { simulationOptMethodController.updatelightForce(LightIntensity.value); });
                LightIntensity.onValueChanged.AddListener(delegate { theLamp.setLightColor(LightIntensity.value, lightMaxIntensityExplicit); });
                LightAngle.onValueChanged.AddListener(delegate { simulationOptMethodController.updatelightAngle(LightAngle.value); });
                LightAngle.onValueChanged.AddListener(delegate { theLamp.setViewingAngle(LightAngle.value); });
                break;
            default:
                break;
        }

        initSliders();
    }

    void initSliders()
    {
        // init slider values.
        switch (simulationMode)
        {
            case SimulationMethod.ExplicitGPU:
                SliderDistance.GetComponent<UpdateValue>().updateSliderBar(CPU3D.nodeDistance, 0.2f, 6f);
                SliderStiffness.GetComponent<UpdateValue>().updateSliderBar(CPU3D.stiffness, 2.0f, 100f);
                SliderGravity.GetComponent<UpdateValue>().updateSliderBar(CPU3D.gravity, 0.05f, 5.0f);
                SliderMaxTravel.GetComponent<UpdateValue>().updateSliderBar(CPU3D.maxTravelDistance, 0.5f, 20f);
                SliderBendingStiffness.GetComponent<UpdateValue>().updateSliderBar(CPU3D.bendingStiffness, 0.01f, 0.5f);
                SliderFriction.GetComponent<UpdateValue>().updateSliderBar(CPU3D.velocityDecay, 0.998f, 0.9999f);

                LightIntensity.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightForce, 0.0f, lightMaxIntensityImplict);
                LightAngle.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightAngle, 1.0f, 45.0f);
                break;
            case SimulationMethod.LocalGlobal:
                SliderDistance.GetComponent<UpdateValue>().updateSliderBar(OptMethod.nodeDistance, 0.2f, 10f);
                SliderStiffness.GetComponent<UpdateValue>().updateSliderBar(OptMethod.stiffness, 2.0f, 100f);
                SliderGravity.GetComponent<UpdateValue>().updateSliderBar(OptMethod.gravity, 0.05f, 10.0f);
                SliderMaxTravel.gameObject.SetActive(false);
                SliderBendingStiffness.gameObject.SetActive(false);
                SliderFriction.gameObject.SetActive(false);
                LightIntensity.GetComponent<UpdateValue>().updateSliderBar(OptMethod.lightForce, 0.0f, lightMaxIntensityExplicit);
                LightAngle.GetComponent<UpdateValue>().updateSliderBar(OptMethod.lightAngle, 1.0f, 45.0f);
                break;
            default:
                break;
        }
    }
}
