using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.UIElements;

// [ExecuteInEditMode]
public class WorldStatus : MonoBehaviour
{
    public SceneNode TheRoot;

    private void Start()
    {

    }

    private void Update()
    {
        Matrix4x4 i = Matrix4x4.identity;
        TheRoot.CompositeXform(ref i);
    }
}
