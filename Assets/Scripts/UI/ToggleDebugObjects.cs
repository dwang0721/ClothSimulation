using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ToggleDebugObjects : MonoBehaviour
{
    public GameObject gameObjectToToggle;
    public Toggle selecteToggle;

    // Start is called before the first frame update
    void Start()
    {
        selecteToggle = gameObject.GetComponent<Toggle>();
        selecteToggle.onValueChanged.AddListener(delegate { toggleGameObjectVisibility(selecteToggle); });

        Debug.Assert(gameObjectToToggle);
        Debug.Assert(selecteToggle);
    }

    void toggleGameObjectVisibility(Toggle tgValue)
    {
        if (tgValue.isOn)
        {
            gameObjectToToggle.SetActive(true);
        }
        else {
            gameObjectToToggle.SetActive(false);
        }
    }

}
