using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class UpdateToggle : MonoBehaviour
{
    Toggle checkBox;
    public int index = -1;
    public UIControlHandler ui;

    // Start is called before the first frame update
    void Start()
    {
        checkBox = gameObject.GetComponent<Toggle>();
        checkBox.onValueChanged.AddListener(delegate { ChangeUICurrentMode(); });
    }

    // Update is called once per frame
    void Update()
    {
        Debug.Assert(checkBox);
        Debug.Assert(index!=-1);
        Debug.Assert(ui);
    }

    void ChangeUICurrentMode()
    {
        if (checkBox.isOn) {
            ui.setCurrentMode(index);
        }        
    }
}
