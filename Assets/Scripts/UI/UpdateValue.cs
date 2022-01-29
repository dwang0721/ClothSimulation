using System.Collections;
using System.Collections.Generic;
using System.Transactions;
using UnityEngine;
using UnityEngine.UI;

// [ExecuteInEditMode]
public class UpdateValue : MonoBehaviour
{
    Text value;
    Slider slider;

    private void Awake()
    {
        slider = gameObject.GetComponent<Slider>();
    }

    // Start is called before the first frame update
    void Start()
    {
        value = gameObject.transform.Find("attrValue").GetComponent<Text>();
    }

    // Update is called once per frame
    void Update()
    {
        value.text = slider.value.ToString("0.00");
    }

    public void updateSliderBar(float curValue, float minValue, float maxValue)
    {
        slider.maxValue = maxValue;
        slider.minValue = minValue;
        slider.value = curValue;
    }
}
