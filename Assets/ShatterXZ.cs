/* Part of the algorithms used in this script are distributed under the MIT License,
 * see below, Copyright (c) 2020 Erik Nordeus. In many of the algorithms, tweaks have
 * been made to better fit the needs of this project.
 * 
 * MIT License
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
using System.Linq;

using LinearAlgebra;
using Voronoi;

namespace VoronoiXZ
{
    [RequireComponent(typeof(MeshFilter))]
    [RequireComponent(typeof(MeshRenderer))]
    [RequireComponent(typeof(MeshCollider))]
    public class ShatterXZ : MonoBehaviour
    {
        public bool showCells = false;
        public bool showTriangles = false;
        public bool showSites = false;
        public bool showVertices = false;

        public bool useSeed = true;
        public int seed = 0;

        public MeshRenderer m_renderer;
        public MeshCollider m_collider;
        public Material m_material;

        public int numberOfPoints = 3;

        private List<Vector3> randomSites = new List<Vector3>();
        public VoronoiDiagram diagram;

        // Boundaries
        private float minX, maxX, minZ, maxZ;
        private float yPos;
        private bool shouldUpdateY = true;
        public void Start()
        {
            m_collider.sharedMesh = GetComponent<MeshFilter>().mesh;
            minX = m_collider.bounds.min.x;
            maxX = m_collider.bounds.max.x;
            minZ = m_collider.bounds.min.z;
            maxZ = m_collider.bounds.max.z;
            yPos = m_collider.bounds.max.y;
        }

        public void Update()
        {
            minX = m_collider.bounds.min.x;
            maxX = m_collider.bounds.max.x;
            minZ = m_collider.bounds.min.z;
            maxZ = m_collider.bounds.max.z;
            yPos = m_collider.bounds.max.y;

            // Always update y-position in case plane is moving
            if (shouldUpdateY)
                yPos = m_collider.bounds.max.y;

            // Show/Hide gizmos (only visible in editor-view)
            if (Input.GetKeyDown(KeyCode.Q))
            {
                if (showSites == true)
                    showSites = false;
                else
                    showSites = true;
            }
            if (Input.GetKeyDown(KeyCode.W))
            {
                if (showTriangles == true)
                    showTriangles = false;
                else
                    showTriangles = true;
            }
            if (Input.GetKeyDown(KeyCode.E))
            {
                if (showVertices == true)
                    showVertices = false;
                else
                    showVertices = true;
            }
            if (Input.GetKeyDown(KeyCode.R))
            {
                if (showCells == true)
                    showCells = false;
                else
                    showCells = true;
            }

            // Manually generate a diagram
            if (Input.GetKeyUp(KeyCode.Space))
            {
                // Generate the voronoi diagram
                GenerateRandomSites(minX, maxX, minZ, maxZ);
                List<VoronoiCell> cells = DelaunayToVoronoi.GenerateVoronoiDiagram(randomSites, minX, maxX, minZ, maxZ, yPos);
                diagram = new VoronoiDiagram(cells);
            }
            // Manually shatter the mesh
            if (Input.GetKeyUp(KeyCode.Return))
            {
                Shatter();
            }
        }

        public void OnCollisionEnter(Collision collision)
        {
            // Stop updating Y so that all pieces are generated at the same height
            shouldUpdateY = false;

            // Generate the voronoi diagram
            GenerateRandomSites(minX, maxX, minZ, maxZ);
            List<VoronoiCell> cells = DelaunayToVoronoi.GenerateVoronoiDiagram(randomSites, minX, maxX, minZ, maxZ, yPos);
            diagram = new VoronoiDiagram(cells);

            Shatter();
        }

        /// <summary>
        /// Generates random sites within the given boundaries. (Boundaries are set in the editor.)
        /// </summary>
        /// <param name="minX">Minimum x-coordinate for boundary.</param>
        /// <param name="maxX">Maximum x-coordinate for boundary.</param>
        /// <param name="minZ">Minimum y-coordinate for boundary.</param>
        /// <param name="maxZ">Maximum y-coordinate for boundary.</param>
        public void GenerateRandomSites(float minX, float maxX, float minZ, float maxZ)
        {

            // Delete any previous sites
            randomSites.Clear();

            if (useSeed)
                Random.InitState(seed);

            // Generate random points
            for (int i = 0; i < numberOfPoints; i++)
            {
                float randomZ = Random.Range(minZ, maxZ);
                float randomX = Random.Range(minX, maxX);
                randomSites.Add(new Vector3(randomX, yPos, randomZ));
            }

            // This makes sure that the star's corners cannot be on the same side
            if (minX > 0)
                minX *= -1;
            if (minZ > 0)
                minZ *= -1;
            if (maxX < 0)
                maxX *= -1;
            if (maxZ < 0)
                maxZ *= -1;

            // Create a big "star" around the random points to handle edge cells.
            randomSites.Add(new Vector3(0f, yPos, minZ * 5f));
            randomSites.Add(new Vector3(0f, yPos, maxZ * 5f));
            randomSites.Add(new Vector3(minX * 5f, yPos, 0f));
            randomSites.Add(new Vector3(maxX * 5f, yPos, 0f));

        }

        /// <summary>
        /// Given a set of points in a cartesian plane, returns the set of points that form a convex hull for the given set.
        /// Based on the Jarvis March algorithm by R. A. Jarvis.
        /// </summary>
        /// <param name="points">Set of points in a cartesian plane</param>
        /// <returns></returns>
        private List<Vertex> CreateConvexHull(List<Vertex> points)
        {
            // Degenerate case
            if (points.Count < 2)
                Debug.Log("At least three (3) points are required for the construction of a convex polygon.");

            // Set of vertices which will form the convex hull
            List<Vertex> vertices = new List<Vertex>();

            // Find the leftmost point in the array which is guaranteed to be on the hull
            Vector3 pointOnHull = points[0].position;
            foreach (Vertex point in points)
            {
                if (point.position.x < pointOnHull.x)
                    pointOnHull = point.position;
            }

            int i = 0;
            while (true)
            {
                vertices.Add(new Vertex(pointOnHull));
                Vector3 endpoint = points[points.Count - 1].position;

                for (int j = 0; j < points.Count; j++)
                {
                    if (endpoint == pointOnHull || IsLeftOfLine(vertices[i].position, endpoint, points[j].position))
                        endpoint = points[j].position;   // Found greater left turn, update
                }

                i++;
                pointOnHull = endpoint;

                // Safety
                if (i > 10000)
                {
                    Debug.Log("Fail at cell");
                    break;
                }

                // End condition
                if (endpoint == vertices[0].position) // Wrapped around to first hull point
                    break;
            }

            return vertices;
        }

        /// <summary>
        /// Returns true if a given point is located to the left of a straight line.
        /// </summary>
        /// <param name="start">Starting point of the line.</param>
        /// <param name="end">End point of the line.</param>
        /// <param name="point">Point of interest.</param>
        /// <returns></returns>
        private bool IsLeftOfLine(Vector3 start, Vector3 end, Vector3 point)
        {
            // d is the signed shortest distance from the line to the point
            float d = ((point.x - start.x) * (end.z - start.z)) - ((point.z - start.z) * (end.x - start.x));

            // Evaluating the distance, if d < 0, then the point is on the left of the line
            if (d < 0)
                return true;
            else
                return false;
        }

        /// <summary>
        /// According to voronoi diagram, creates new game objects in the shape of cells.
        /// </summary>
        private void Shatter()
        {
            for (int i = 0; i < diagram.cells.Count; i++)
            {
                List<Vertex> vertices = diagram.cells[i].GetVertices();

                // These are the helper cells, we don't want to render those.
                if (diagram.cells[i].sitePosition == randomSites[randomSites.Count - 4] || diagram.cells[i].sitePosition == randomSites[randomSites.Count - 3] || diagram.cells[i].sitePosition == randomSites[randomSites.Count - 2] || diagram.cells[i].sitePosition == randomSites[randomSites.Count - 1])
                {
                    Debug.Log("Removed helper cell at site: " + diagram.cells[i].sitePosition);
                    continue;
                }
                // If there, for some reason, are less than 3 vertices in our cell, we cannot create a mesh from it.
                else if (vertices.Count < 3)
                {
                    Debug.Log("Less than three vertices in cell: " + i);
                    continue;
                }

                // Sort the vertices by creating a convex hull. This is required for the clipping algorithm.
                List<Vertex> convexHull = CreateConvexHull(vertices);

                // Define our cut polygon (the plane).
                List<Vector3> clipPlane = new List<Vector3>();
                clipPlane.Add(new Vector3(minX, yPos, minZ));
                clipPlane.Add(new Vector3(maxX, yPos, minZ));
                clipPlane.Add(new Vector3(maxX, yPos, maxZ));
                clipPlane.Add(new Vector3(minX, yPos, maxZ));

                // Clip cell
                List<Vertex> polygon = Intersections.ClipPolygon(convexHull, clipPlane);

                // Triangulate the new polygon
                List<Triangle> triangles = IncrementalTriangulationAlgorithm.TriangulatePoints(polygon);
                for (int k = 0; k < triangles.Count; k++)
                {
                    // Make sure all faces are facing the same way
                    if (!DelaunayTriangulationAlgorithm.IsTriangleOrientedClockwise(triangles[k].v1.GetPos2D_XZ(), triangles[k].v2.GetPos2D_XZ(), triangles[k].v3.GetPos2D_XZ()))
                    {
                        triangles[k].ChangeOrientation();
                    }
                }

                Mesh piece_mesh = TriangleMesh(triangles);
                GameObject piece = new GameObject("cell" + i);
                piece.AddComponent<MeshFilter>();
                piece.AddComponent<MeshRenderer>();
                piece.AddComponent<Rigidbody>();

                piece.GetComponent<MeshRenderer>().material = m_material;
                piece.GetComponent<MeshFilter>().mesh = piece_mesh;

                piece.AddComponent<MeshCollider>();
                piece.GetComponent<MeshCollider>().convex = true;
                piece.GetComponent<MeshCollider>().sharedMesh = piece_mesh;
            }
            Destroy(this.gameObject);
        }

        /// <summary>
        /// Creates a triangle mesh from a given list of triangles. Guarantees that
        /// each vertex exists only once in the final mesh. Inspiration from 
        /// https://docs.unity3d.com/ScriptReference/Mesh.html
        /// </summary>
        /// <param name="triangles"></param>
        /// <returns>A triangle Mesh from the given list of triangles</returns>
        public static Mesh TriangleMesh(List<Triangle> triangles)
        {
            // If no exist, than non exists. Just move on 
            if (triangles == null)
            {
                return null;
            }

            // Create the list with unique vertices. This, and all that follows, will be of Vector3 in order to set the vertices of the mesh
            //A list will make it fast to check if a vertex already exists in the collection
            List<Vector3> newVertices = new List<Vector3>();

            foreach (Triangle t in triangles)
            {
                Vector3 v1 = t.v1.position;
                Vector3 v2 = t.v2.position;
                Vector3 v3 = t.v3.position;

                if (!newVertices.Contains(v1))
                {
                    newVertices.Add(v1);
                }
                if (!newVertices.Contains(v2))
                {
                    newVertices.Add(v2);
                }
                if (!newVertices.Contains(v3))
                {
                    newVertices.Add(v3);
                }
            }

            //Create the list to add all vertices for the mesh
            List<Vector3> meshVertices = new List<Vector3>(newVertices);

            //Create the list to add triangles by using the vertices for the mesh
            List<int> meshTriangles = new List<int>();

            //Use a dictionay to quickly find which positon in the list a Vector3 has
            Dictionary<Vector3, int> verticePos = new Dictionary<Vector3, int>();

            for (int i = 0; i < meshVertices.Count; i++)
            {
                verticePos.Add(meshVertices[i], i);
            }

            foreach (Triangle t in triangles)
            {
                Vector3 v1 = t.v1.position;
                Vector3 v2 = t.v2.position;
                Vector3 v3 = t.v3.position;

                meshTriangles.Add(verticePos[v1]);
                meshTriangles.Add(verticePos[v2]);
                meshTriangles.Add(verticePos[v3]);
            }


            // Create the final mesh
            Mesh mesh = new Mesh();
            Vector3[] meshVerticesArray = new Vector3[meshVertices.Count];

            for (int i = 0; i < meshVerticesArray.Length; i++)
            {
                Vector3 v = meshVertices[i];

                meshVerticesArray[i] = new Vector3(v.x, v.y, v.z);
            }

            mesh.vertices = meshVerticesArray;
            mesh.triangles = meshTriangles.ToArray();

            mesh.RecalculateNormals();

            return mesh;
        }

        /// <summary>
        /// Displays the voronoi sites using gizmos spheres. (Only visible in the editor view).
        /// </summary>
        /// <param name="p">List of sites to display.</param>
        /// <param name="radius">Raidus of draw spheres.</param>
        /// <param name="color">Color of draw spheres.</param>
        public static void DisplaySites(List<Vector3> p, float radius, Color color)
        {

            if (p == null) return;

            Gizmos.color = color;

            foreach (Vector3 points in p)
            {
                Gizmos.DrawSphere(points, radius);
            }
        }

        /// <summary>
        /// Draw the triangle edges in the mesh using gizmos draw lines. (Only visible in editor view).
        /// </summary>
        /// <param name="mesh">Mesh to draw.</param>
        /// <param name="color">Color of draw lines.</param>
        public static void DisplayEdge(Mesh mesh, Color color)
        {
            //If no mesh exist, just move on. 
            if (mesh == null)
            {
                return;
            }

            // The tringles and vertices from the mesh
            int[] meshTriangles = mesh.triangles;
            Vector3[] meshVertices = mesh.vertices;

            for (int i = 0; i < meshTriangles.Length; i += 3)
            {
                // Get the mesh vertices from the triangles and move three steps for three vertices needed
                Vector3 v1 = meshVertices[meshTriangles[i + 0]];
                Vector3 v2 = meshVertices[meshTriangles[i + 1]];
                Vector3 v3 = meshVertices[meshTriangles[i + 2]];

                Gizmos.color = color;

                //Display the triangles and use Gizmos.Drawline. The line is between the vertices.
                Gizmos.DrawLine(v1, v2);
                Gizmos.DrawLine(v2, v3);
                Gizmos.DrawLine(v3, v1);
            }
        }

        public static void DisplayVoronoiEdge(List<VoronoiEdge> edges, Color color)
        {
            if (edges.Count <= 0)
                return;

            foreach (VoronoiEdge edge in edges)
            {
                Gizmos.color = color;
                Gizmos.DrawLine(edge.v1, edge.v2);
            }
        }

        /// <summary>
        /// Draw the cells of a voronoi diagram with a random colour using gizmos.
        /// </summary>
        /// <param name="cells">List of cells in the voronoi diagrams.</param>
        public void DisplayVoronoiCells(List<VoronoiCell> cells)
        {
            Random.InitState(seed);

            for (int i = 0; i < cells.Count; i++)
            {
                VoronoiCell c = cells[i];
                Vector3 p1 = c.sitePosition;

                Gizmos.color = new Color(Random.Range(0f, 1f), Random.Range(0f, 1f), Random.Range(0f, 1f));

                List<Vector3> vertices = new List<Vector3>();
                List<int> triangles = new List<int>();

                vertices.Add(p1);
                for (int j = 0; j < c.edges.Count; j++)
                {
                    Vector3 p3 = c.edges[j].v1;
                    Vector3 p2 = c.edges[j].v2;

                    vertices.Add(p2);
                    vertices.Add(p3);

                    triangles.Add(0);
                    triangles.Add(vertices.Count - 2);
                    triangles.Add(vertices.Count - 1);
                }

                Mesh triangleMesh = new Mesh();

                triangleMesh.vertices = vertices.ToArray();
                triangleMesh.triangles = triangles.ToArray();

                triangleMesh.RecalculateNormals();

                Gizmos.DrawMesh(triangleMesh);
            }
        }

        /// <summary>
        /// Draws the different components of the voronoi diagram using gizmos. (Only visible in the editor view.)
        /// </summary>
        private void OnDrawGizmos()
        {
            // Display voronoi cells
            if (showCells)
            {
                DisplayVoronoiCells(diagram.cells);

            }

            if (showTriangles)
            {
                // Display triangles
                Mesh meshDelaunay = TriangleMesh(DelaunayTriangulationAlgorithm.DelaunayTriangulation(randomSites));
                DisplayEdge(meshDelaunay, Color.blue);
            }

            if (showSites)
                // Dislpay sites
                DisplaySites(randomSites, 0.1f, Color.blue);

            if (showVertices)
            {

                // Display voronoi edges
                for (int i = 0; i < diagram.cells.Count; i++)
                {
                    DisplayVoronoiEdge(diagram.cells[i].GetEdges(), Color.red);
                }
                List<Vector3> voronoiVertices = new List<Vector3>();
                foreach (VoronoiCell cell in diagram.cells)
                {
                    List<Vertex> temp = cell.GetVertices();
                    foreach (Vertex v in temp)
                        voronoiVertices.Add(v.position);
                }

                DisplaySites(voronoiVertices, 0.1f, Color.red);
            }

            // Display boundaries
            Gizmos.color = Color.black;
            Gizmos.DrawLine(new Vector3(maxX, yPos, minZ), new Vector3(minX, yPos, minZ));
            Gizmos.DrawLine(new Vector3(maxX, yPos, maxZ), new Vector3(maxX, yPos, minZ));
            Gizmos.DrawLine(new Vector3(minX, yPos, maxZ), new Vector3(maxX, yPos, maxZ));
            Gizmos.DrawLine(new Vector3(minX, yPos, minZ), new Vector3(minX, yPos, maxZ));
        }

    }

    public static class IncrementalTriangulationAlgorithm
    {
        /// <summary>
        /// Triangulates a given set of points using an incremental triangulation algorithm.
        /// </summary>
        /// <param name="points">Set of points to triangulate.</param>
        /// <returns>List of triangles generated from the set of points.</returns>
        public static List<Triangle> TriangulatePoints(List<Vertex> points)
        {
            List<Triangle> triangles = new List<Triangle>();

            // Sort the points along x-axis in increasing order.
            points = points.OrderBy(n => n.position.x).ToList();

            // The first 3 vertices always form a triangle
            Triangle newTriangle = new Triangle(points[0].position, points[1].position, points[2].position);
            triangles.Add(newTriangle);

            // All edges that form the triangles
            List<Edge> edges = new List<Edge>();

            edges.Add(new Edge(newTriangle.v1, newTriangle.v2));
            edges.Add(new Edge(newTriangle.v2, newTriangle.v3));
            edges.Add(new Edge(newTriangle.v3, newTriangle.v1));

            // Now add the remaining triangles one by one
            for (int i = 3; i < points.Count; i++)
            {
                Vector3 currentPoint = points[i].position;
                List<Edge> newEdges = new List<Edge>();

                // Is the edge visible? We only need to check if the midpoint of the edge is visible
                for (int j = 0; j < edges.Count; j++)
                {
                    Edge currentEdge = edges[j];
                    Vector3 midPoint = (currentEdge.v1.position + currentEdge.v2.position) / 2f;

                    Edge edgeToMidpoint = new Edge(currentPoint, midPoint);

                    bool canSeeEdge = true;
                    for (int k = 0; k < edges.Count; k++)
                    {
                        // Don't compare edge with itself
                        if (k == j)
                            continue;
                        // If edges are intersecting, don't
                        if (AreEdgesIntersecting(edgeToMidpoint, edges[k]))
                        {
                            canSeeEdge = false;
                            break;
                        }
                    }

                    // If we've made it this far, we have a valid triangle
                    if (canSeeEdge)
                    {
                        Edge edgeToPoint1 = new Edge(currentEdge.v1, new Vertex(currentPoint));
                        Edge edgeToPoint2 = new Edge(currentEdge.v2, new Vertex(currentPoint));

                        newEdges.Add(edgeToPoint1);
                        newEdges.Add(edgeToPoint2);

                        newTriangle = new Triangle(edgeToPoint1.v1, edgeToPoint1.v2, edgeToPoint2.v1);

                        triangles.Add(newTriangle);
                    }
                }

                for (int j = 0; j < newEdges.Count; j++)
                {
                    edges.Add(newEdges[j]);
                }
            }

            return triangles;
        }

        /// <summary>
        /// Returns true if two given edges are intersecting.
        /// </summary>
        /// <param name="edge1">First edge.</param>
        /// <param name="edge2">Second edge.</param>
        /// <returns>Returns wether the edges are intersecting or not.</returns>
        private static bool AreEdgesIntersecting(Edge edge1, Edge edge2)
        {
            Vector2 l1_p1 = new Vector2(edge1.v1.position.x, edge1.v1.position.z);
            Vector2 l1_p2 = new Vector2(edge1.v2.position.x, edge1.v2.position.z);

            Vector2 l2_p1 = new Vector2(edge2.v1.position.x, edge2.v1.position.z);
            Vector2 l2_p2 = new Vector2(edge2.v2.position.x, edge2.v2.position.z);

            bool isIntersecting = Intersections.AreLinesIntersecting(l1_p1, l1_p2, l2_p1, l2_p2);

            return isIntersecting;
        }
    }

    /// <summary>
    /// Class holding functions needed for a delaunay triangulation algorithm.
    /// </summary>
    public class DelaunayTriangulationAlgorithm
    {
        /// <summary>
        /// Converts a list of triangles into half edges.
        /// </summary>
        /// <param name="triangles">List of triangles to convert.</param>
        /// <returns>Returns a list of half edges generated from given triangles.</returns>
        public static List<HalfEdge> Triangles2HalfEdges(List<Triangle> triangles)
        {
            // Make sure all triangles have the same orientation (clockwise
            foreach (Triangle triangle in triangles)
            {
                Vector2 p1 = triangle.v1.GetPos2D_XZ();
                Vector2 p2 = triangle.v2.GetPos2D_XZ();
                Vector2 p3 = triangle.v3.GetPos2D_XZ();

                if (!IsTriangleOrientedClockwise(p1, p2, p3))
                    triangle.ChangeOrientation();
            }

            // Construct the half edges from all the triangles
            List<HalfEdge> halfEdges = new List<HalfEdge>(triangles.Count * 3);
            foreach (Triangle triangle in triangles)
            {
                HalfEdge he1 = new HalfEdge(triangle.v1);
                HalfEdge he2 = new HalfEdge(triangle.v2);
                HalfEdge he3 = new HalfEdge(triangle.v3);

                he1.nextEdge = he2;
                he2.nextEdge = he3;
                he3.nextEdge = he1;

                he1.prevEdge = he3;
                he2.prevEdge = he1;
                he3.prevEdge = he2;


                he1.v.halfEdge = he2;
                he2.v.halfEdge = he3;
                he3.v.halfEdge = he1;

                triangle.halfEdge = he1;

                he1.t = triangle;
                he2.t = triangle;
                he3.t = triangle;

                halfEdges.Add(he1);
                halfEdges.Add(he2);
                halfEdges.Add(he3);
            }

            // Find the opposite half edge for all half edges
            for (int i = 0; i < halfEdges.Count; i++)
            {
                HalfEdge he = halfEdges[i];

                Vertex goingToVertex = he.v;
                Vertex goingFromVertex = he.prevEdge.v;

                for (int j = 0; j < halfEdges.Count; j++)
                {
                    if (i == j)
                        continue;

                    HalfEdge heOpposite = halfEdges[j];

                    // Break if we find the opposite half edge of current half edge
                    if (goingFromVertex.position == heOpposite.v.position && goingToVertex.position == heOpposite.prevEdge.v.position)
                    {
                        he.oppositeEdge = heOpposite;
                        break;
                    }
                }
            }

            return halfEdges;
        }

        /// <summary>
        /// Returns true if a triangle is oriented clockwise.
        /// </summary>
        /// <param name="p1">First vertex of triangle.</param>
        /// <param name="p2">Second vertex of triangle.</param>
        /// <param name="p3">Third vertex of triangle.</param>
        /// <returns>Returns wether a triangle is oriented clockwise or not.</returns>
        public static bool IsTriangleOrientedClockwise(Vector2 p1, Vector2 p2, Vector2 p3)
        {
            bool isClockwise = true;
            float determinant = p1.x * p2.y + p3.x * p1.y + p2.x * p3.y - p1.x * p3.y - p3.x * p2.y - p2.x * p1.y;

            // If determinant is greater than zero, then it is not clockwise and should change orientation
            if (determinant > 0f)
                isClockwise = false;

            return isClockwise;
        }

        /// <summary>
        /// Check where a point lies in relation to a circle.
        /// </summary>
        /// <param name="aVec">First point on the circle.</param>
        /// <param name="bVec">Second point on the circle.</param>
        /// <param name="cVec">Third point on the circle.</param>
        /// <param name="dVec">Point to compare to circle.</param>
        /// <returns></returns>
        public static float PointInRelationToCircle(Vector2 aVec, Vector2 bVec, Vector2 cVec, Vector2 dVec)
        {
            float a = aVec.x - dVec.x;
            float d = bVec.x - dVec.x;
            float g = cVec.x - dVec.x;

            float b = aVec.y - dVec.y;
            float e = bVec.y - dVec.y;
            float h = cVec.y - dVec.y;

            float c = a * a + b * b;
            float f = d * d + e * e;
            float i = g * g + h * h;

            float determinant = (a * e * i) + (b * f * g) + (c * d * h) - (g * e * c) - (h * f * a) - (i * d * b);

            return determinant;
        }

        /// <summary>
        /// Returns true if a quadrilateral is convex. (Assumes no 3 points are colinear and no hourglass shapes.)
        /// </summary>
        /// <param name="a">First vertex of quadrilateral.</param>
        /// <param name="b">Second vertex of quadrilateral.</param>
        /// <param name="c">Third vertex of quadrilateral.</param>
        /// <param name="d">Fourth vertex of quadrilateral.</param>
        /// <returns>Returns wether or not the quadrilateral is convex.</returns>
        public static bool IsQaudrilateralConvex(Vector2 a, Vector2 b, Vector2 c, Vector2 d)
        {
            bool isConvex = false;

            bool abc = IsTriangleOrientedClockwise(a, b, c);
            bool abd = IsTriangleOrientedClockwise(a, b, d);
            bool bcd = IsTriangleOrientedClockwise(b, c, d);
            bool cad = IsTriangleOrientedClockwise(c, a, d);

            if (abc && abd && bcd & !cad)
                isConvex = true;
            else if (abc && abd && !bcd & cad)
                isConvex = true;
            else if (abc && !abd && bcd & cad)
                isConvex = true;
            else if (!abc && !abd && !bcd & cad)
                isConvex = true;
            else if (!abc && !abd && bcd & !cad)
                isConvex = true;
            else if (!abc && abd && !bcd & !cad)
                isConvex = true;

            return isConvex;
        }

        /// <summary>
        /// Generates a delaunay triangulation from a given set of points.
        /// </summary>
        /// <param name="sites">Set of points to triangulate.</param>
        /// <returns>A list of triangles generated from the given set of points.</returns>
        public static List<Triangle> DelaunayTriangulation(List<Vector3> sites)
        {
            // 1. Triangulate given set of points with arbitrary triangulation algorithm
            List<Vertex> vertices = new List<Vertex>();
            for (int i = 0; i < sites.Count; i++)
            {
                vertices.Add(new Vertex(sites[i]));
            }

            List<Triangle> triangles = IncrementalTriangulationAlgorithm.TriangulatePoints(vertices);

            // 2. In order to speed up edge flipping, convert triangles to half edges since those are easier to flip.
            List<HalfEdge> halfEdges = Triangles2HalfEdges(triangles);

            // 3. Flip edges until we have delaunay triangulation
            int safety = 0; // Keep track of iterations
            int flippedEdges = 0;   // Keep track of how many edges have been flipped

            while (true)
            {
                safety += 1;

                if (safety > 100000)
                {
                    Debug.Log("Stuck in an endless loop");
                    break;
                }

                bool hasFlippedEdge = false;

                // Search through all edges to see if we can flip an edge
                for (int i = 0; i < halfEdges.Count; i++)
                {
                    HalfEdge thisEdge = halfEdges[i];

                    // Is this edge sharing an edge, otherwise it is a border, and then we cannot flip
                    if (thisEdge.oppositeEdge == null)
                    {
                        continue;
                    }

                    Vertex a = thisEdge.v;
                    Vertex b = thisEdge.nextEdge.v;
                    Vertex c = thisEdge.prevEdge.v;
                    Vertex d = thisEdge.oppositeEdge.nextEdge.v;

                    Vector2 aPos = a.GetPos2D_XZ();
                    Vector2 bPos = b.GetPos2D_XZ();
                    Vector2 cPos = c.GetPos2D_XZ();
                    Vector2 dPos = d.GetPos2D_XZ();

                    // Use the circle test to test if we need to flip this edge
                    if (PointInRelationToCircle(aPos, bPos, cPos, dPos) < 0f)
                    {
                        // Are these the two triangles that share this edge forming a convex quadrilateral
                        // If not, edges cannot be flipped
                        if (IsQaudrilateralConvex(aPos, bPos, cPos, dPos))
                        {
                            // If the new triangle after a flip is not better, then don't flip
                            // This will also stop the algorithm from ending up in an endless loop
                            if (PointInRelationToCircle(bPos, cPos, dPos, aPos) < 0f)
                            {
                                continue;
                            }

                            flippedEdges += 1;
                            hasFlippedEdge = true;
                            FlipEdge(thisEdge);
                        }
                    }
                }

                // We have searched through all edges and have not found an edge to flip, we have found a valid triangulation
                if (!hasFlippedEdge)
                {
                    break;
                }
            }
            return triangles;
        }

        /// <summary>
        /// Flips a given half edge. See https://en.wikipedia.org/wiki/Delaunay_triangulation for visual representation.
        /// Note: This function assumes that the half edge meets the flipping conditions, use only if those conditions are met.
        /// </summary>
        /// <param name="one">The first half edge of the flip.</param>
        public static void FlipEdge(HalfEdge one)
        {
            /*** Required data ***/
            // Together with one, these edges form one's triangle
            HalfEdge two = one.nextEdge;
            HalfEdge three = one.prevEdge;

            // The opposite half edge's triangle
            HalfEdge four = one.oppositeEdge;
            HalfEdge five = one.oppositeEdge.nextEdge;
            HalfEdge six = one.oppositeEdge.prevEdge;

            // Vertices
            Vertex a = one.v;
            Vertex b = one.nextEdge.v;
            Vertex c = one.prevEdge.v;
            Vertex d = one.oppositeEdge.nextEdge.v;


            /*** Flip ***/

            a.halfEdge = one.nextEdge;
            c.halfEdge = one.oppositeEdge.nextEdge;

            one.nextEdge = three;
            one.prevEdge = five;

            two.nextEdge = four;
            two.prevEdge = six;

            three.nextEdge = five;
            three.prevEdge = one;

            four.nextEdge = six;
            four.prevEdge = two;

            five.nextEdge = one;
            five.prevEdge = three;

            six.nextEdge = two;
            six.prevEdge = four;

            one.v = b;
            two.v = b;
            three.v = c;
            four.v = d;
            five.v = d;
            six.v = a;

            Triangle t1 = one.t;
            Triangle t2 = four.t;

            one.t = t1;
            three.t = t1;
            five.t = t1;

            two.t = t2;
            four.t = t2;
            six.t = t2;

            t1.v1 = b;
            t1.v2 = c;
            t1.v3 = d;

            t2.v1 = b;
            t2.v2 = d;
            t2.v3 = a;

            t1.halfEdge = three;
            t2.halfEdge = four;

        }
    }
    /// <summary>
    /// Class holding necessary functions to generate a voronoi diagram from a delaunay triangulation.
    /// </summary>
    public class DelaunayToVoronoi
    {
        /// <summary>
        /// Generates a voronoi diagram from a given set of points using delaunay triangulation.
        /// </summary>
        /// <param name="sites">Set of points.</param>
        /// <returns>Returns a list containing the voronoi cells of the generated voronoi diagram.</returns>
        public static List<VoronoiCell> GenerateVoronoiDiagram(List<Vector3> sites, float minX, float maxX, float minZ, float maxZ, float yPos)
        {
            List<Triangle> triangles = DelaunayTriangulationAlgorithm.DelaunayTriangulation(sites);

            List<VoronoiEdge> voronoiEdges = new List<VoronoiEdge>();

            for (int i = 0; i < triangles.Count; i++)
            {
                Triangle t = triangles[i];

                HalfEdge e1 = t.halfEdge;
                HalfEdge e2 = e1.nextEdge;
                HalfEdge e3 = e2.nextEdge;

                Vector3 v1 = e1.v.position;
                Vector3 v2 = e2.v.position;
                Vector3 v3 = e3.v.position;

                Vector2 v12D = new Vector2(v1.x, v1.z);
                Vector2 v22D = new Vector2(v2.x, v2.z);
                Vector2 v32D = new Vector2(v3.x, v3.z);

                Vector2 center2D = CalculateCircleCenter(v12D, v22D, v32D);

                Vector3 voronoiVertex = new Vector3(center2D.x, yPos, center2D.y);

                TryAddingVoronoiEdgeFromTriangleEdge(e1, voronoiVertex, voronoiEdges, minX, maxX, minZ, maxZ);
                TryAddingVoronoiEdgeFromTriangleEdge(e2, voronoiVertex, voronoiEdges, minX, maxX, minZ, maxZ);
                TryAddingVoronoiEdgeFromTriangleEdge(e3, voronoiVertex, voronoiEdges, minX, maxX, minZ, maxZ);
            }

            List<VoronoiCell> voronoiCells = new List<VoronoiCell>();

            for (int i = 0; i < voronoiEdges.Count; i++)
            {
                VoronoiEdge e = voronoiEdges[i];

                int cellPos = TryFindCellPos(e, voronoiCells);

                if (cellPos == -1)
                {
                    VoronoiCell newCell = new VoronoiCell(e.sitePosition);
                    voronoiCells.Add(newCell);
                    newCell.edges.Add(e);
                }
                else
                {
                    voronoiCells[cellPos].edges.Add(e);
                }
            }

            return voronoiCells;
        }

        /// <summary>
        /// Calculate the center that intersects all three given points. Note: currently only supports 2D.
        /// </summary>
        /// <param name="p1">First point.</param>
        /// <param name="p2">Second point.</param>
        /// <param name="p3">Third point.</param>
        /// <returns>Returns a 2D vector containing the coordinates for the circle's center.</returns>
        private static Vector2 CalculateCircleCenter(Vector2 p1, Vector2 p2, Vector2 p3)
        {
            //Debug.Log(p1 + ", " + p2 + ", " + p3);
            Vector2 center = new Vector2();

            float ma = (p2.y - p1.y) / (p2.x - p1.x);
            float mb = (p3.y - p2.y) / (p3.x - p2.x);

            center.x = (ma * mb * (p1.y - p3.y) + mb * (p1.x + p2.x) - ma * (p2.x + p3.x)) / (2 * (mb - ma));
            center.y = (-1 / ma) * (center.x - (p1.x + p2.x) / 2) + (p1.y + p2.y) / 2;

            return center;
        }

        /// <summary>
        /// Checks if there already exists a voronoi site, and uses it if it does. If not, returns -1.
        /// </summary>
        /// <param name="e">Voronoi edge to check against.</param>
        /// <param name="voronoiCells">List of voronoi cells in the voronoi diagram.</param>
        /// <returns>Returns the index of existing site if possible, else -1.</returns>
        private static int TryFindCellPos(VoronoiEdge e, List<VoronoiCell> voronoiCells)
        {
            for (int i = 0; i < voronoiCells.Count; i++)
            {
                if (e.sitePosition == voronoiCells[i].sitePosition)
                    return i;
            }

            return -1;
        }
        /// <summary>
        /// Tries to add a voronoi edge from a triangle edge if possible. (Not possible if there is no opposite half edge.)
        /// </summary>
        /// <param name="e">Half edge of triangle.</param>
        /// <param name="voronoiVertex">Voronoi vertex.</param>
        /// <param name="allEdges">List containting all voronoi edges in the voronoi diagram.</param>
        private static void TryAddingVoronoiEdgeFromTriangleEdge(HalfEdge e, Vector3 voronoiVertex, List<VoronoiEdge> allEdges, float minX, float maxX, float minZ, float maxZ)
        {
            Vector2 center2D;
            // If there is no triangle on the other side, no edge can be created, connect it to the closest corner instead.
            if (e.oppositeEdge == null)
            {
                return;
                //center2D = FindClosestCorner(voronoiVertex, minX, maxX, minZ, maxZ);
            }
            else
            {
                // However, if there is a triangle on the other side, we can connect its circumcenter with our edge.
                HalfEdge eNeighbour = e.oppositeEdge;

                Vector3 v1 = eNeighbour.v.position;
                Vector3 v2 = eNeighbour.nextEdge.v.position;
                Vector3 v3 = eNeighbour.nextEdge.nextEdge.v.position;

                Vector2 v12D = new Vector2(v1.x, v1.z);
                Vector2 v22D = new Vector2(v2.x, v2.z);
                Vector2 v32D = new Vector2(v3.x, v3.z);

                center2D = CalculateCircleCenter(v12D, v22D, v32D);
            }

            Vector3 voronoiVertexNeighbour = new Vector3(center2D.x, voronoiVertex.y, center2D.y);

            // Create the new edge and add it to the edge list.
            VoronoiEdge edge = new VoronoiEdge(voronoiVertex, voronoiVertexNeighbour, e.prevEdge.v.position);
            if (edge.v1 == edge.v2)
                return;

            allEdges.Add(edge);
        }
    }

    /// <summary>
    /// Class holding functions for handling intersections.
    /// </summary>
    public class Intersections
    {
        /// <summary>
        /// Check if two lines are intersecting.
        /// </summary>
        /// <param name="l1_p1">Line 1 point 1.</param>
        /// <param name="l1_p2">Line 1 point 2.</param>
        /// <param name="l2_p1">Line 2 point 1.</param>
        /// <param name="l2_p2">Line 2 point 2.</param>
        /// <returns>Returns true if the lines are intersecting, else returns false.</returns>
        public static bool AreLinesIntersecting(Vector2 l1_p1, Vector2 l1_p2, Vector2 l2_p1, Vector2 l2_p2)
        {
            bool isIntersecting = false;

            float den = (l2_p2.y - l2_p1.y) * (l1_p2.x - l1_p1.x) - (l2_p2.x - l2_p1.x) * (l1_p2.y - l1_p1.y);

            if (den != 0f)
            {
                float u_a = ((l2_p2.x - l2_p1.x) * (l1_p1.y - l2_p1.y) - (l2_p2.y - l2_p1.y) * (l1_p1.x - l2_p1.x)) / den;
                float u_b = ((l1_p2.x - l1_p1.x) * (l1_p1.y - l2_p1.y) - (l1_p2.y - l1_p1.y) * (l1_p1.x - l2_p1.x)) / den;

                if (u_a >= 0f && u_a <= 1f && u_b >= 0f && u_b <= 1f)
                    isIntersecting = true;
            }

            return isIntersecting;
        }

        public static void IntersectionPoint(Vector2 l1_p1, Vector2 l1_p2, Vector2 l2_p1, Vector2 l2_p2, List<Vector2> intersections)
        {
            // Get direction of lines
            Vector2 l1_dir = (l1_p2 - l1_p1).normalized;
            Vector2 l2_dir = (l2_p2 - l2_p1).normalized;

            // Get normals of lines
            Vector2 l1_normal = new Vector2(-l1_dir.y, l1_dir.x);
            Vector2 l2_normal = new Vector2(-l2_dir.y, l2_dir.x);

            // Get normal form of lines (Ax + By = k1, Cx + Dy = k2)
            float A = l1_normal.x;
            float B = l1_normal.y;

            float C = l2_normal.x;
            float D = l2_normal.y;

            float k1 = (A * l1_p1.x) + (B * l1_p1.y);
            float k2 = (A * l2_p1.x) + (B * l2_p1.y);

            // Check degenerate cases
            if (IsParallel(l1_normal, l2_normal))
                // Lines are parallel, no solutions
                return;

            if (IsOrthogonal(l1_p1 - l2_p1, l1_normal))
                // Lines are the same, infinite solutions
                return;

            // Calculate intersection
            float x_intersect = (D * k1 - B * k2) / (A * D - B * C);
            float y_intersect = (-C * k1 + A * k2) / (A * D - B * C);
            Vector2 intersectionPoint = new Vector2(x_intersect, y_intersect);

            // Check if point is within line segments
            if (IsBetween(l1_p1, l1_p2, intersectionPoint) && IsBetween(l2_p1, l2_p2, intersectionPoint))
            {
                Debug.Log("Intersecting edge: " + l1_p1 + ", " + l1_p2);
                Debug.Log("Intersection point: " + intersectionPoint.x + ", " + intersectionPoint.y);
                intersections.Add(intersectionPoint);
            }

        }

        /// <summary>
        /// Checks wether two vectors are parallel or not.
        /// </summary>
        /// <param name="v1">First vector.</param>
        /// <param name="v2">Second vector.</param>
        /// <returns>True if parallel, false if not.</returns>
        public static bool IsParallel(Vector2 v1, Vector2 v2)
        {
            //2 vectors are parallel if the angle between the vectors are 0 or 180 degrees
            if (Vector2.Angle(v1, v2) == 0f || Vector2.Angle(v1, v2) == 180f)
            {
                return true;
            }

            return false;
        }
        /// <summary>
        /// Checks wether two vectors are orthogonal two each other.
        /// </summary>
        /// <param name="v1">First vector.</param>
        /// <param name="v2">Second vector.</param>
        /// <returns>True if orthogonal, false if not.</returns>
        public static bool IsOrthogonal(Vector2 v1, Vector2 v2)
        {
            //2 vectors are orthogonal is the dot product is 0
            //We have to check if close to 0 because of floating numbers
            if (Mathf.Abs(Vector2.Dot(v1, v2)) < 0.000001f)
            {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Checks wether a given point is between two points on a line.
        /// </summary>
        /// <param name="a">First point on line.</param>
        /// <param name="b">Second point on line.</param>
        /// <param name="c">Point of interest.</param>
        /// <returns>True if between, false it not.</returns>
        public static bool IsBetween(Vector2 a, Vector2 b, Vector2 c)
        {
            bool isBetween = false;

            //Entire line segment
            Vector2 ab = b - a;
            //The intersection and the first point
            Vector2 ac = c - a;

            //Need to check 2 things: 
            //1. If the vectors are pointing in the same direction = if the dot product is positive
            //2. If the length of the vector between the intersection and the first point is smaller than the entire line
            if (Vector2.Dot(ab, ac) > 0f && ab.sqrMagnitude >= ac.sqrMagnitude)
            {
                isBetween = true;
            }

            return isBetween;
        }

        /// <summary>
        /// Returns the vertices defining the intersection of polygon 1 and polygon 2 according to the Sutherland-Hodgman algorithm.
        /// </summary>
        /// <param name="poly_1">The polygon to cut.</param>
        /// <param name="poly_2">The stencil polygon.</param>
        /// <returns>List of vertices defining the intersection of the polygons.</returns>
        public static List<Vertex> ClipPolygon(List<Vertex> poly_1, List<Vector3> poly_2)
        {
            List<LinearAlgebra.Plane> clippingPlanes = new List<LinearAlgebra.Plane>();

            for (int i = 0; i < poly_2.Count; i++)
            {
                int iPlusOne = ((i + 1 % poly_2.Count) + poly_2.Count) % poly_2.Count;

                Vector3 v1 = poly_2[i];
                Vector3 v2 = poly_2[iPlusOne];

                Vector3 planePos = (v1 + v2) * 0.5f;

                Vector3 planeDir = v2 - v1;

                // Should point inwards
                Vector3 planeNormal = new Vector3(-planeDir.z, 0f, planeDir.x).normalized;

                clippingPlanes.Add(new LinearAlgebra.Plane(planePos, planeNormal));

            }

            List<Vertex> vertices = ClipPolygon(poly_1, clippingPlanes);

            return vertices;
        }

        /// <summary>
        /// Cuts a polygon according to clipping planes according to the Sutherland-Hodgman algorithm.
        /// </summary>
        /// <param name="poly_1">Polygon to cut.</param>
        /// <param name="clippingPlanes">List of clipping planes to use.</param>
        /// <returns>List of vertices defining the new polyon.</returns>
        public static List<Vertex> ClipPolygon(List<Vertex> poly_1, List<LinearAlgebra.Plane> clippingPlanes)
        {
            // Clone the vertices since the algorithm will remove some
            List<Vertex> vertices = new List<Vertex>(poly_1);

            // Save the new vertices temporarily before transfering them to final list
            List<Vertex> temp_vertices = new List<Vertex>();

            // Clipping tajm
            foreach (LinearAlgebra.Plane plane in clippingPlanes)
            {
                for (int i = 0; i < vertices.Count; i++)
                {
                    int iPlusOne = ((i + 1 % vertices.Count) + vertices.Count) % vertices.Count;

                    Vertex v1 = vertices[i];
                    Vertex v2 = vertices[iPlusOne];

                    // Calculate the distance to the plane for each vertex, since the planes will point inwards, a positive distance means inside.
                    float distance_v1 = Vector3.Dot(plane.normal, v1.position - plane.pos);
                    float distance_v2 = Vector3.Dot(plane.normal, v2.position - plane.pos);
                    // If both points are outside, do nothing
                    // If both points are inside, keep v2
                    if (distance_v1 > 0f && distance_v2 > 0f)
                        temp_vertices.Add(v2);
                    // v1 outside, v2 inside, save intersection point and v2
                    else if (distance_v1 < 0f && distance_v2 > 0f)
                    {
                        Vector3 rayDir = (v2.position - v1.position).normalized;
                        Vector3 intersectionPoint = GetRayPlaneIntersection(plane.pos, plane.normal, v1.position, rayDir);

                        temp_vertices.Add(new Vertex(intersectionPoint));
                        temp_vertices.Add(v2);
                    }
                    // v1 inside, v2 outside, save intersection point
                    else if (distance_v1 > 0f && distance_v2 < 0f)
                    {
                        Vector3 rayDir = (v2.position - v1.position).normalized;
                        Vector3 intersectionPoint = GetRayPlaneIntersection(plane.pos, plane.normal, v1.position, rayDir);

                        temp_vertices.Add(new Vertex(intersectionPoint));
                    }
                }

                vertices.Clear();
                vertices.AddRange(temp_vertices);
                temp_vertices.Clear();
            }

            return vertices;
        }

        /// <summary>
        /// Returns the point of intersection between a ray and a plane.
        /// </summary>
        /// <param name="planePos">Point in plane.</param>
        /// <param name="planeNormal">Normal vector of plane.</param>
        /// <param name="rayStart">Starting position of ray.</param>
        /// <param name="rayDir">Ending position of ray.</param>
        /// <returns>Intersection point as a Vector3</returns>
        public static Vector3 GetRayPlaneIntersection(Vector3 planePos, Vector3 planeNormal, Vector3 rayStart, Vector3 rayDir)
        {
            float den = Vector3.Dot(-planeNormal, rayDir);

            Vector3 vecBetween = planePos - rayStart;

            float t = Vector3.Dot(vecBetween, -planeNormal) / den;

            Vector3 intersectionPoint = rayStart + rayDir * t;

            return intersectionPoint;
        }
    }
}
