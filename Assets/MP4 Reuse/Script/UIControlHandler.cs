using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class UIControlHandler: MonoBehaviour
{
    // reference to the UI
    public enum ControlMode
    {
        Translation,
        Rotation,
        Scale
    }

    // reference to the Scene
    public GameObject theLamp;
    GameObject theHead;
    public CPU3D simulationController3D;
    ControlMode currentMode;

    // reference to all the UI elements
    Toggle checkBoxT, checkBoxR, checkBoxS;
    Slider sliderX, sliderY, sliderZ;

    private void Awake()
    {
    }

    // Start is called before the first frame update
    void Start()
    {
        theHead = theLamp.GetComponent<Lamp>().head.gameObject;
        // Initialization.
        checkBoxT = gameObject.transform.Find("CheckBoxT").GetComponent<Toggle>();
        checkBoxR = gameObject.transform.Find("CheckBoxR").GetComponent<Toggle>();
        checkBoxS = gameObject.transform.Find("CheckBoxS").GetComponent<Toggle>();

        sliderX = gameObject.transform.Find("SliderX").GetComponent<Slider>();
        sliderY = gameObject.transform.Find("SliderY").GetComponent<Slider>();
        sliderZ = gameObject.transform.Find("SliderZ").GetComponent<Slider>();
        sliderX.onValueChanged.AddListener(delegate { updateSelectedItem(); });
        sliderY.onValueChanged.AddListener(delegate { updateSelectedItem(); });
        sliderZ.onValueChanged.AddListener(delegate { updateSelectedItem(); });
        setCurrentMode(0);
    }

    // Update is called once per frame
    void Update()
    {
        Debug.Assert(checkBoxT & checkBoxR & checkBoxS);
        Debug.Assert(sliderX & sliderY & sliderZ);
        Debug.Assert(theLamp);
        Debug.Assert(theHead);
    }

    /// UI Controls
    public void updateSelectedItem()
    {
        switch (currentMode)
        {
            case ControlMode.Translation:
                theLamp.transform.localPosition = new Vector3(sliderX.value, sliderY.value, sliderZ.value);
                break;
            case ControlMode.Rotation:
                theHead.transform.localEulerAngles = new Vector3(sliderX.value, sliderY.value, 0.0f);
                break;
            case ControlMode.Scale:
                break;
        }
        //simulationController3D.updateLightDirection(theHead.transform.position, );
    }

    public void setCurrentMode(int index)
    {
        switch (index)
        {
            case 0:
                currentMode = ControlMode.Translation;
                sliderX.gameObject.SetActive(true);
                sliderY.gameObject.SetActive(true);
                sliderZ.gameObject.SetActive(true);
                break;
            case 1:
                currentMode = ControlMode.Rotation;
                sliderX.gameObject.SetActive(true);
                sliderY.gameObject.SetActive(true);
                sliderZ.gameObject.SetActive(false);
                break;
            case 2:
                currentMode = ControlMode.Scale;
                sliderX.gameObject.SetActive(false);
                sliderY.gameObject.SetActive(false);
                sliderZ.gameObject.SetActive(false);
                break;
        }
        updateSliders();
    }

    public void updateSliders()
    {
        Debug.Assert(sliderX & sliderY & sliderZ);
        Vector3 values = new Vector3 (0.0f, 0.0f, 0.0f);
 
        // find max and min value
        switch (currentMode)
        {
            case ControlMode.Translation:
                values.x = theLamp.transform.localPosition.x;
                values.y = theLamp.transform.localPosition.y;
                values.z = theLamp.transform.localPosition.z;
                sliderX.GetComponent<UpdateValue>().updateSliderBar(values.x, -40, 40);
                sliderY.GetComponent<UpdateValue>().updateSliderBar(values.y, -40, 40);
                sliderZ.GetComponent<UpdateValue>().updateSliderBar(values.z, -60, -10);
                break;
            case ControlMode.Rotation:
                float angleX = theHead.transform.localRotation.eulerAngles.x;
                float angleY = theHead.transform.localRotation.eulerAngles.y;

                angleX = angleX > 180.0f ? -(360f - angleX) : angleX;
                angleX = angleX < -180.0f ? (360f + angleX) : angleX;
                angleY = angleY > 180.0f ? -(360f - angleY) : angleY;
                angleY = angleY < -180.0f ? (360f + angleY) : angleY;

                sliderX.GetComponent<UpdateValue>().updateSliderBar(angleX, -80, 80);
                sliderY.GetComponent<UpdateValue>().updateSliderBar(angleY, -80, 80);
                break;
            case ControlMode.Scale:
                break;
        }        
    }
}
