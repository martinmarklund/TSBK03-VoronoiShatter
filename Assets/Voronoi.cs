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


using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using LinearAlgebra;

namespace Voronoi
{
    [System.Serializable]
    public class VoronoiDiagram
    {
        public List<VoronoiCell> cells;

        public VoronoiDiagram(List<VoronoiCell> cells)
        {
            this.cells = cells;
        }

        public List<VoronoiCell> GetCells()
        {
            return cells;
        }
    }

    /// <summary>
    /// Class for voronoi edge, contains information about vertices and which site if belongs to.
    /// </summary>
    [System.Serializable]
    public class VoronoiEdge
    {
        public Vector3 v1;
        public Vector3 v2;

        public Vector3 sitePosition;

        public VoronoiEdge(Vector3 v1, Vector3 v2, Vector3 sitePosition)
        {
            this.v1 = v1;
            this.v2 = v2;

            this.sitePosition = sitePosition;
        }
    }

    /// <summary>
    /// Class for voronoi cell, contains information about site position and its edges.
    /// </summary>
    [System.Serializable]
    public class VoronoiCell
    {
        public Vector3 sitePosition;

        public List<VoronoiEdge> edges = new List<VoronoiEdge>();

        public VoronoiCell(Vector3 sitePosition)
        {
            this.sitePosition = sitePosition;
        }

        /// <summary>
        /// Get the edges of the cell.
        /// </summary>
        /// <returns>Returns a list of the cells edges.</returns>
        public List<VoronoiEdge> GetEdges()
        {
            return edges;
        }

        public List<Vertex> GetVertices()
        {
            List<Vertex> vertices = new List<Vertex>();
            for (int i = 0; i < edges.Count; i++)
            {
                vertices.Add(new Vertex(edges[i].v1));
            }

            return vertices;
        }
    }
}
