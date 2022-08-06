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
    public SimulationScenario simulationScenario;
    Slider SliderStiffness, SliderGravity, SliderDistance, SliderMaxTravel, SliderBendingStiffness, SliderFriction;
    Slider SliderLightIntensity, SliderLightAngle, SliderExFriction;

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
        simulationScenario = globalData.simulationScenario;

        Debug.Log(string.Format("Simualtion Configuration: Mode: {0}, Scenario: {1}, Resolution: {2}", simulationMode, simulationScenario, globalData.resolution));

        // set up the run time simulator
        switch (simulationMode)
        {
            case SimulationMethod.ExplicitGPU:
                simulationController3D.gameObject.SetActive(true);
                simulationOptMethodController.gameObject.SetActive(false);

                simulationController3D.resolution = globalData.resolution;
                simulationController3D.simulationScenario = globalData.simulationScenario;
                simulationController3D.clothDebugNodeSize = globalData.clothDebugNodeSize;
                simulationController3D.nodeDistance = globalData.nodeDistance;
                simulationController3D.collisionPushAwayDistance = globalData.collisionPushAwayDistance;
                simulationController3D.gravity = globalData.gravity;
                simulationController3D.stiffness = globalData.stiffness;
                simulationController3D.maxTravelDistance = globalData.maxTravelDistance;
                simulationController3D.bendingStiffness = globalData.bendingStiffness;
                simulationController3D.velocityDecay = globalData.velocityDecay;

                theLamp.gameObject.SetActive(true);
                theLamp.simulationMode = SimulationMethod.ExplicitGPU;
                break;
            case SimulationMethod.LocalGlobal:
                simulationController3D.gameObject.SetActive(false);
                simulationOptMethodController.gameObject.SetActive(true);

                simulationOptMethodController.resolution = globalData.resolution;
                simulationOptMethodController.simulationScenario = globalData.simulationScenario;
                simulationOptMethodController.clothDebugNodeSize = globalData.clothDebugNodeSize;
                simulationOptMethodController.nodeDistance = globalData.nodeDistance * 2;
                simulationOptMethodController.collisionPushAwayDistance = globalData.collisionPushAwayDistance;
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

        SliderLightIntensity = gameObject.transform.Find("RunTimeConfigPanel").Find("intensity").GetComponent<Slider>();
        SliderLightAngle = gameObject.transform.Find("RunTimeConfigPanel").Find("angle").GetComponent<Slider>();
        SliderExFriction = gameObject.transform.Find("RunTimeConfigPanel").Find("extFriction").GetComponent<Slider>();

        Debug.Assert(SliderDistance && SliderStiffness && SliderGravity && SliderMaxTravel && SliderBendingStiffness && SliderFriction);
        Debug.Assert(SliderLightIntensity && SliderLightAngle);

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

                SliderLightIntensity.onValueChanged.AddListener(delegate { simulationController3D.updatelightForce(SliderLightIntensity.value); });
                SliderLightIntensity.onValueChanged.AddListener(delegate { theLamp.setLightColor(SliderLightIntensity.value, lightMaxIntensityImplict); });
                SliderLightAngle.onValueChanged.AddListener(delegate { simulationController3D.updatelightAngle(SliderLightAngle.value); });
                SliderLightAngle.onValueChanged.AddListener(delegate { theLamp.setViewingAngle(SliderLightAngle.value); });
                SliderExFriction.onValueChanged.AddListener(delegate { simulationController3D.setExtFriction(SliderExFriction.value); });
                break;
            case SimulationMethod.LocalGlobal:
                SliderDistance.onValueChanged.AddListener(delegate { simulationOptMethodController.updateNodeDistance(SliderDistance.value); });
                SliderStiffness.onValueChanged.AddListener(delegate { simulationOptMethodController.updateStiffness(SliderStiffness.value); });
                SliderGravity.onValueChanged.AddListener(delegate { simulationOptMethodController.updateGravity(SliderGravity.value); });

                SliderLightIntensity.onValueChanged.AddListener(delegate { simulationOptMethodController.updatelightForce(SliderLightIntensity.value); });
                SliderLightIntensity.onValueChanged.AddListener(delegate { theLamp.setLightColor(SliderLightIntensity.value, lightMaxIntensityExplicit); });
                SliderLightAngle.onValueChanged.AddListener(delegate { simulationOptMethodController.updatelightAngle(SliderLightAngle.value); });
                SliderLightAngle.onValueChanged.AddListener(delegate { theLamp.setViewingAngle(SliderLightAngle.value); });
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
                SliderDistance.GetComponent<UpdateValue>().updateSliderBar(globalData.nodeDistance, 0.2f, 6f);
                SliderStiffness.GetComponent<UpdateValue>().updateSliderBar(globalData.stiffness, 2.0f, 100f);
                SliderGravity.GetComponent<UpdateValue>().updateSliderBar(globalData.gravity, 0.05f, 5.0f);
                SliderMaxTravel.GetComponent<UpdateValue>().updateSliderBar(globalData.maxTravelDistance, 0.5f, 20f);
                SliderBendingStiffness.GetComponent<UpdateValue>().updateSliderBar(globalData.bendingStiffness, 0.01f, 0.5f);
                SliderFriction.GetComponent<UpdateValue>().updateSliderBar(globalData.velocityDecay, 0.998f, 0.9999f);

                SliderLightIntensity.GetComponent<UpdateValue>().updateSliderBar(globalData.lightForce, 0.0f, lightMaxIntensityImplict);
                SliderLightAngle.GetComponent<UpdateValue>().updateSliderBar(globalData.lightAngle, 1.0f, 45.0f);
                SliderExFriction.GetComponent<UpdateValue>().updateSliderBar(globalData.extFriction, 0.0f, 1.0f);
                break;
            case SimulationMethod.LocalGlobal:
                SliderDistance.GetComponent<UpdateValue>().updateSliderBar(globalData.nodeDistance, 0.2f, 10f);
                SliderStiffness.GetComponent<UpdateValue>().updateSliderBar(globalData.stiffness, 2.0f, 100f);
                SliderGravity.GetComponent<UpdateValue>().updateSliderBar(globalData.gravity, 0.05f, 10.0f);
                SliderMaxTravel.gameObject.SetActive(false);
                SliderBendingStiffness.gameObject.SetActive(false);
                SliderFriction.gameObject.SetActive(false);
                SliderLightIntensity.GetComponent<UpdateValue>().updateSliderBar(globalData.lightForce, 0.0f, lightMaxIntensityExplicit);
                SliderLightAngle.GetComponent<UpdateValue>().updateSliderBar(globalData.lightAngle, 1.0f, 45.0f);
                SliderExFriction.gameObject.SetActive(false);
                break;
            default:
                break;
        }
    }
}
