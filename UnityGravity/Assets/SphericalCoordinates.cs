using UnityEngine;


public class SphericalCoordinates
{
    public float Radius = 10.0f;
    public float Inclination = Mathf.PI / 2;
    public float Rotation = 0.0f;

    public SphericalCoordinates() { }
    public SphericalCoordinates(float radius, float inclination, float rotation)
    {
        Radius = radius;
        Inclination = inclination;
        Rotation = rotation;
    }

    public Vector3 toCartesian
    {
        get
        {
            float a = Radius * Mathf.Cos(Inclination);
            return new Vector3(a * Mathf.Cos(Rotation), Radius * Mathf.Sin(Inclination), a * Mathf.Sin(Rotation));
        }
    }
}