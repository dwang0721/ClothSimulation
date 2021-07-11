using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GizmoNode : MonoBehaviour
{

    protected Matrix4x4 mCombinedParentXform;
    public SceneNode mSelectedSceneNode;
    public Vector3 NodeOrigin = Vector3.zero;
    public List<NodePrimitive> PrimitiveList;

    // Use this for initialization
    protected void Start()
    {
        InitializeSceneNode();
    }

    // Update is called once per frame
    void Update()
    {
        mCombinedParentXform = mSelectedSceneNode.getCombinedParentXform();
        CompositeXform(ref mCombinedParentXform);
    }

    public void SetReferenceSceneNode(SceneNode selectedSceneNode)
    {
        mSelectedSceneNode = selectedSceneNode;
    }

    private void InitializeSceneNode()
    {
        SetReferenceSceneNode(mSelectedSceneNode);
        mCombinedParentXform = mSelectedSceneNode.getCombinedParentXform();
    }

    // This must be called _BEFORE_ each draw!! 
    public void CompositeXform(ref Matrix4x4 parentXform)
    {
        Matrix4x4 orgT = Matrix4x4.Translate(NodeOrigin);
        Matrix4x4 trs = Matrix4x4.TRS(transform.localPosition, transform.localRotation, transform.localScale);

        mCombinedParentXform = parentXform * orgT * trs;

        // propagate to all children
        foreach (Transform child in transform)
        {
            SceneNode cn = child.GetComponent<SceneNode>();
            if (cn != null)
            {
                cn.CompositeXform(ref mCombinedParentXform);
            }
        }

        // disenminate to primitives
        foreach (NodePrimitive p in PrimitiveList)
        {
            p.LoadShaderMatrix(ref mCombinedParentXform);
        }
    }
}
