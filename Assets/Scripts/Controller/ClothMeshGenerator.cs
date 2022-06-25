using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

public class ClothMeshGenerator : MonoBehaviour
{
    Mesh clothMesh;
    private int clothResolutionX, clothResolutionY;
    private void Start()
    {
        clothMesh = GetComponent<MeshFilter>().mesh;
    }

    void buildClothMesh()
    {
        clothMesh.Clear();

        int nCol = clothResolutionX;
        int nRow = clothResolutionY;

        int nVetices = nCol * nRow;
        int nTriangle = 2 * (nCol - 1) * (nRow - 1);
        int nQuads = 2 * (nCol - 1) * (nRow - 1);

        Vector3[] v = new Vector3[nVetices];
        Vector3[] n = new Vector3[nVetices];
        Vector2[] uvs = new Vector2[nVetices];
        int[] t = new int[nQuads * 3];

        int triIndex = 0;
        int verIndex = 0;

        float s = 6.0f + Mathf.Sin(Time.time); ;

        v[0] = new Vector3(-s, s, 0);
        v[1] = new Vector3(s, s, 0);
        v[2] = new Vector3(-s, -s, 0);
        v[3] = new Vector3(s, -s, 0);

        n[0] = new Vector3(0, 0, -1);
        n[1] = new Vector3(0, 0, -1);
        n[2] = new Vector3(0, 0, -1);
        n[3] = new Vector3(0, 0, -1);

        uvs[0] = new Vector2(0, 1);
        uvs[1] = new Vector2(1, 1);
        uvs[2] = new Vector2(0, 0);
        uvs[3] = new Vector2(1, 0);

        t[triIndex] = verIndex;
        t[triIndex + 1] = verIndex + nCol + 1;
        t[triIndex + 2] = verIndex + nCol;
        t[triIndex + 3] = verIndex;
        t[triIndex + 4] = verIndex + 1;
        t[triIndex + 5] = verIndex + nCol + 1;

        clothMesh.vertices = v;
        clothMesh.triangles = t;
        clothMesh.normals = n;
        clothMesh.uv = uvs;
    }

    // initMesh set up pos, uv, normal, and index once for the initial setup
    public void initMesh(Vector3[] vertices, Vector3[] normals, Vector2[] uvs, int[] triangles)
    {
        clothMesh.vertices = vertices;
        clothMesh.normals = normals;
        clothMesh.triangles = triangles;
        clothMesh.uv = uvs;
    }

    public void initMeshExplict(HairNode3D[] clothNodesArray, int dim_x, int dim_y)
    {
        // dim_x is number of hairs 
        // dim_y is number of nodes per hair
        // clothNodesArray iterating in y direction. 
        
        clothResolutionX = dim_x;
        clothResolutionY = dim_y;

        Debug.Assert(clothNodesArray.Length == dim_x * dim_y);
        Debug.Assert(clothMesh);

        int nVertices = dim_x * dim_y;
        int nTriangles = 2 * (dim_x - 1) * (dim_y - 1);
        int nTriangleIndices = 3 * nTriangles;

        Vector3[] vert = new Vector3[nVertices];
        Vector3[] norm = new Vector3[nVertices];
        Vector2[] uvs = new Vector2[nVertices];
        int[] triangle = new int[nTriangleIndices];
        int triIndex = 0;

        for (int i = 0; i < nVertices; i++)
        {
            vert[i] = new Vector3(clothNodesArray[i].x, clothNodesArray[i].y, clothNodesArray[i].z);
            norm[i] = new Vector3(0, 0, 1);

            int x = i % dim_x;
            int y = i / dim_x;
            float u = x * 1.0f / (dim_x - 1);
            float v = 1.0f - y * 1.0f / (dim_y - 1);
            uvs[i] = new Vector2(u, v);

            if (x < (dim_x - 1) && y < (dim_y - 1))
            {
                triangle[triIndex] = i;
                triangle[triIndex + 1] = i + dim_x;
                triangle[triIndex + 2] = i + dim_x + 1;
                triangle[triIndex + 3] = i;
                triangle[triIndex + 4] = i + dim_x + 1;
                triangle[triIndex + 5] = i + 1;
                triIndex += 6;
            }
        }

        clothMesh.vertices = vert;
        clothMesh.triangles = triangle;
        clothMesh.normals = norm;
        clothMesh.uv = uvs;
    }

    public void initMeshImplicit(Matrix<double> position, int dim_x, int dim_y)
    {
        Debug.Assert(position.RowCount == 3* dim_x * dim_y);
        Debug.Assert(clothMesh);

        clothResolutionX = dim_x;
        clothResolutionY = dim_y;

        int nVertices = dim_x * dim_y;
        int nTriangles = 2 * (dim_x - 1) * (dim_y - 1);
        int nTriangleIndices = 3 * nTriangles;

        Vector3[] vert = new Vector3[nVertices];
        Vector3[] norm = new Vector3[nVertices];
        Vector2[] uvs = new Vector2[nVertices];
        int[] triangle = new int[nTriangleIndices];
        int triIndex = 0;

        for (int i = 0; i < nVertices; i++)
        {
            vert[i] = new Vector3((float)(position.At(3 * i, 0)),
                                  (float)(position.At(3 * i + 1, 0)),
                                  (float)(position.At(3 * i + 2, 0)));
            norm[i] = new Vector3(0, 0, 1);

            int x = i % dim_x;
            int y = i / dim_x;
            float u = x * 1.0f / (dim_x - 1);
            float v = 1.0f - y * 1.0f / (dim_y - 1);
            uvs[i] = new Vector2(u, v);

            if (x < (dim_x - 1) && y < (dim_y - 1))
            {
                triangle[triIndex] = i;
                triangle[triIndex + 1] = i + dim_x;
                triangle[triIndex + 2] = i + dim_x + 1;
                triangle[triIndex + 3] = i;
                triangle[triIndex + 4] = i + dim_x + 1;
                triangle[triIndex + 5] = i + 1;
                triIndex += 6;
            }

            clothMesh.vertices = vert;
            clothMesh.triangles = triangle;
            clothMesh.normals = norm;
            clothMesh.uv = uvs;
        }
    }

        // updateMesh only update the mesh's pos and normal. 
        public void updateMeshExplict(HairNode3D[] clothNodesArray)
    {
        //clothMesh.Clear();

        Debug.Assert(clothResolutionX * clothResolutionY == clothNodesArray.Length);

        int nVertices = clothResolutionX * clothResolutionY;
        Vector3[] vert = new Vector3[nVertices];
        Vector3[] norm = new Vector3[nVertices];
        for (int i = 0; i < nVertices; i++)
        {
            vert[i] = new Vector3(clothNodesArray[i].x, clothNodesArray[i].y, clothNodesArray[i].z);
            norm[i] = new Vector3(clothNodesArray[i].nx, clothNodesArray[i].ny, clothNodesArray[i].nz);
        }

        clothMesh.vertices = vert;
        clothMesh.normals = norm;
    }

    public void updateMeshImplict(Matrix<double> position)
    {
        //clothMesh.Clear();

        Debug.Assert(clothResolutionX * clothResolutionY * 3 == position.RowCount);

        int nVertices = clothResolutionX * clothResolutionY;
        Vector3[] vert = new Vector3[nVertices];
        for (int i = 0; i < nVertices; i++)
        {
            vert[i] = new Vector3((float)(position.At(3 * i, 0)),
                                  (float)(position.At(3 * i + 1, 0)),
                                  (float)(position.At(3 * i + 2, 0)));
        }

        clothMesh.vertices = vert;
    }
}
