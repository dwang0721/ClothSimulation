using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public enum SimulationMethod { ExplicitGPU, ExplicitCPU, LocalGlobal };

public class SimulationControl : MonoBehaviour
{
    public CPU3D simulationController3D;
    public OptMethod simulationOptMethodController;
    public Lamp theLamp;
    public LampOptMethod theLampOptMethod;
    Slider SliderStiffness, SliderGravity, SliderDistance, SliderMaxTravel, SliderBendingStiffness, SliderFriction;
    Slider LightIntensity, LightAngle;

    public SimulationMethod simulationMode = SimulationMethod.LocalGlobal;

    private void Awake()
    {
        switch (simulationMode)
        {
            case SimulationMethod.ExplicitGPU:
                simulationController3D.gameObject.SetActive(true);
                simulationOptMethodController.gameObject.SetActive(false);
                theLamp.gameObject.SetActive(true);
                break;
            case SimulationMethod.LocalGlobal:
                simulationController3D.gameObject.SetActive(false);
                simulationOptMethodController.gameObject.SetActive(true);
                theLamp.gameObject.SetActive(false);
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
        SliderDistance = gameObject.transform.Find("distance").GetComponent<Slider>();
        SliderStiffness = gameObject.transform.Find("stiffness").GetComponent<Slider>();
        SliderGravity = gameObject.transform.Find("gravity").GetComponent<Slider>();        
        SliderMaxTravel = gameObject.transform.Find("maxTravel").GetComponent<Slider>();
        SliderBendingStiffness = gameObject.transform.Find("bendStiff").GetComponent<Slider>();
        SliderFriction = gameObject.transform.Find("friction").GetComponent<Slider>();
        LightIntensity = gameObject.transform.Find("intensity").GetComponent<Slider>();
        LightAngle = gameObject.transform.Find("angle").GetComponent<Slider>();

        // set up callbacks
        if (simulationController3D.gameObject.active)
        {
            SliderDistance.onValueChanged.AddListener(delegate { simulationController3D.updateNodeDistance(SliderDistance.value); });
            SliderStiffness.onValueChanged.AddListener(delegate { simulationController3D.updateStiffness(SliderStiffness.value); });
            SliderGravity.onValueChanged.AddListener(delegate { simulationController3D.updateGravity(SliderGravity.value); });
            SliderMaxTravel.onValueChanged.AddListener(delegate { simulationController3D.updateMaxTravel(SliderMaxTravel.value); });
            SliderBendingStiffness.onValueChanged.AddListener(delegate { simulationController3D.updateBendStiff(SliderBendingStiffness.value); });
            SliderFriction.onValueChanged.AddListener(delegate { simulationController3D.updateFriction(SliderFriction.value); });

            LightIntensity.onValueChanged.AddListener(delegate { simulationController3D.updatelightForce(LightIntensity.value); });
            LightIntensity.onValueChanged.AddListener(delegate { theLamp.setLightColor(LightIntensity.value); });
            LightAngle.onValueChanged.AddListener(delegate { simulationController3D.updatelightAngle(LightAngle.value); });
            LightAngle.onValueChanged.AddListener(delegate { theLamp.setViewingAngle(LightAngle.value); });
        }
        else {
            SliderDistance.onValueChanged.AddListener(delegate { simulationOptMethodController.updateNodeDistance(SliderDistance.value); });
            SliderStiffness.onValueChanged.AddListener(delegate { simulationOptMethodController.updateStiffness(SliderStiffness.value); });
            SliderGravity.onValueChanged.AddListener(delegate { simulationOptMethodController.updateGravity(SliderGravity.value); });

            LightIntensity.onValueChanged.AddListener(delegate { simulationOptMethodController.updatelightForce(LightIntensity.value); });
            LightIntensity.onValueChanged.AddListener(delegate { theLampOptMethod.setLightColor(LightIntensity.value); });
            LightAngle.onValueChanged.AddListener(delegate { simulationOptMethodController.updatelightAngle(LightAngle.value); });
            LightAngle.onValueChanged.AddListener(delegate { theLampOptMethod.setViewingAngle(LightAngle.value); });
        }

        initSliders();
    }

    void initSliders()
    {
        // init slider values.
        if (simulationController3D.gameObject.active)
        {
            SliderDistance.GetComponent<UpdateValue>().updateSliderBar(CPU3D.nodeDistance, 0.2f, 6f);
            SliderStiffness.GetComponent<UpdateValue>().updateSliderBar(CPU3D.stiffness, 2.0f, 100f);
            SliderGravity.GetComponent<UpdateValue>().updateSliderBar(CPU3D.gravity, 0.05f, 5.0f);
            SliderMaxTravel.GetComponent<UpdateValue>().updateSliderBar(CPU3D.maxTravelDistance, 0.5f, 20f);
            SliderBendingStiffness.GetComponent<UpdateValue>().updateSliderBar(CPU3D.bendingStiffness, 0.01f, 0.5f);
            SliderFriction.GetComponent<UpdateValue>().updateSliderBar(CPU3D.velocityDecay, 0.998f, 0.9999f);

            LightIntensity.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightForce, 0.0f, 0.3f);
            LightAngle.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightAngle, 1.0f, 45.0f);
        }
        else {
            SliderDistance.GetComponent<UpdateValue>().updateSliderBar(OptMethod.nodeDistance, 0.2f, 10f);
            SliderStiffness.GetComponent<UpdateValue>().updateSliderBar(OptMethod.stiffness, 2.0f, 100f);
            SliderGravity.GetComponent<UpdateValue>().updateSliderBar(OptMethod.gravity, 0.05f, 10.0f);
            SliderMaxTravel.gameObject.SetActive(false);
            SliderBendingStiffness.gameObject.SetActive(false);
            SliderFriction.gameObject.SetActive(false);

            LightIntensity.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightForce, 0.0f, 0.3f);
            LightAngle.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightAngle, 1.0f, 45.0f);
        }
    }
}
