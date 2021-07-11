using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class SimulationControl : MonoBehaviour
{
    public CPU3D simulationController3D;
    public Lamp theLamp;
    Slider SliderStiffness, SliderGravity, SliderDistance, SliderMaxTravel, SliderBendingStiffness, SliderFriction;
    Slider LightIntensity, LightAngle;
    
    // Start is called before the first frame update
    void Start()
    {
        SliderDistance = gameObject.transform.Find("distance").GetComponent<Slider>();
        SliderStiffness = gameObject.transform.Find("stiffness").GetComponent<Slider>();
        SliderGravity = gameObject.transform.Find("gravity").GetComponent<Slider>();        
        SliderMaxTravel = gameObject.transform.Find("maxTravel").GetComponent<Slider>();
        SliderBendingStiffness = gameObject.transform.Find("bendStiff").GetComponent<Slider>();
        SliderFriction = gameObject.transform.Find("friction").GetComponent<Slider>();

        LightIntensity = gameObject.transform.Find("intensity").GetComponent<Slider>();
        LightAngle = gameObject.transform.Find("angle").GetComponent<Slider>();

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

        initSliders();
    }

    void initSliders()
    {
        SliderDistance.GetComponent<UpdateValue>().updateSliderBar(CPU3D.nodeDistance, 0.2f, 2f);
        SliderStiffness.GetComponent<UpdateValue>().updateSliderBar(CPU3D.stiffness, 3f, 10f);
        SliderGravity.GetComponent<UpdateValue>().updateSliderBar(CPU3D.gravity, 0.05f, 0.3f);
        SliderMaxTravel.GetComponent<UpdateValue>().updateSliderBar(CPU3D.maxTravelDistance, 0.5f, 20f);
        SliderBendingStiffness.GetComponent<UpdateValue>().updateSliderBar(CPU3D.bendingStiffness, 0.01f, 0.5f);
        SliderFriction.GetComponent<UpdateValue>().updateSliderBar(CPU3D.velocityDecay, 0.998f, 0.9999f);

        LightIntensity.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightForce, 0.0f, 0.3f);
        LightAngle.GetComponent<UpdateValue>().updateSliderBar(CPU3D.lightAngle, 1.0f, 45.0f);
    }
}
