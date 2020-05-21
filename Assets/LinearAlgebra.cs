/* MIT License
 *
 * Copyright (c) 2020 Erik Nordeus
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. */

using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace LinearAlgebra
{
    public class Vertex
    {
        public Vector3 position;

        public HalfEdge halfEdge;

        public Triangle triangle;

        public Vertex(Vector3 position)
        {
            this.position = position;
        }

        public Vector2 GetPos2D_XZ()
        {
            Vector2 pos_2d_xz = new Vector2(position.x, position.z);

            return pos_2d_xz;
        }
        
        public Vector2 GetPos2D_XY()
        {
            Vector2 pos_2d_xy = new Vector2(position.x, position.y);

            return pos_2d_xy;
        }

        public Vector2 GetPos2D_ZY()
        {
            Vector2 pos_2d_yz = new Vector2(position.z, position.y);

            return pos_2d_yz;
        }

    }

    public class HalfEdge
    {
        public Vertex v;

        public Triangle t;

        public HalfEdge nextEdge;

        public HalfEdge prevEdge;

        public HalfEdge oppositeEdge;

        public HalfEdge(Vertex v)
        {
            this.v = v;
        }
    }

    public class Triangle
    {
        // Corners
        public Vertex v1;
        public Vertex v2;
        public Vertex v3;

        // If we are using the half edge mesh structure, we just need one half edge
        public HalfEdge halfEdge;

        public Triangle(Vertex v1, Vertex v2, Vertex v3)
        {
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
        }

        public Triangle(Vector3 v1, Vector3 v2, Vector3 v3)
        {
            this.v1 = new Vertex(v1);
            this.v2 = new Vertex(v2);
            this.v3 = new Vertex(v3);
        }

        public Triangle(HalfEdge halfEdge)
        {
            this.halfEdge = halfEdge;
        }

        public void ChangeOrientation()
        {
            Vertex temp = this.v1;

            this.v1 = this.v2;
            this.v2 = temp;
        }
    }

    public class Edge
    {
        public Vertex v1;
        public Vertex v2;

        public bool isIntersecting = false;

        public Edge(Vertex v1, Vertex v2)
        {
            this.v1 = v1;
            this.v2 = v2;
        }

        public Edge(Vector3 v1, Vector3 v2)
        {
            this.v1 = new Vertex(v1);
            this.v2 = new Vertex(v2);
        }

        public Vector2 GetVertex2D(Vertex v)
        {
            return new Vector2(v.position.x, v.position.z);
        }

        public void FlipEdge()
        {
            Vertex temp = v1;
            v1 = v2;
            v2 = temp;
        }
    }

    public class Plane
    {
        public Vector3 pos;
        public Vector3 normal;

        public Plane(Vector3 pos, Vector3 normal)
        {
            this.pos = pos;
            this.normal = normal;
        }
    }
}
