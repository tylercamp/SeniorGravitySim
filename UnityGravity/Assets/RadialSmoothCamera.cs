using UnityEngine;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
//using UnityEditor;
using UnityEngine.Rendering;

public class RadialSmoothCamera : MonoBehaviour
{
    private Vector3? lastPosition = null;
    private List<Vector3> runningMouseDeltas = new List<Vector3>();

    public SphericalCoordinates RadialOrientation = new SphericalCoordinates(10.0f, 0.0f, Mathf.PI / 2);
    

    public float RotateScale = 100.0f;
    public float ScrollScale = 50.0f;
    public float Distance = 5000.0f;

    public Vector3 CenterPosition = Vector3.zero;

    private float rotFriction = 1e-4f;
    private float incFriction = 5e-2f;
    private float incVelocity = 0.0f;
    private float rotVelocity = 0.0f;

	// Use this for initialization
	void Start ()
	{

	}

    void ProcessMouse()
    {
        if (Input.GetMouseButton(0))
        {
            incVelocity *= 0.8f;
            rotVelocity *= 0.8f;

            if (lastPosition.HasValue)
            {
                var delta = Input.mousePosition - lastPosition.Value;
                delta = -delta;
                runningMouseDeltas.Add(delta);
                if (runningMouseDeltas.Count > 60)
                    runningMouseDeltas.RemoveAt(0);

                //if (delta == Vector3.zero)
                //    return;

                RadialOrientation.Inclination += delta.y / (RotateScale * 2);
                RadialOrientation.Rotation += delta.x / (RotateScale);
            }

            lastPosition = Input.mousePosition;
        }
        else
        {
            if (runningMouseDeltas.Count > 0 && runningMouseDeltas.Last() != Vector3.zero)
            {
                var sum = Vector3.zero;
                foreach (var p in runningMouseDeltas)
                    sum += p;
                var avg = sum / runningMouseDeltas.Count;

                rotVelocity = avg.x;
                incVelocity = avg.y;
            }

            lastPosition = null;
            runningMouseDeltas.Clear();
        }

        RadialOrientation.Rotation += rotVelocity / (RotateScale);
        RadialOrientation.Inclination += incVelocity / (RotateScale * 2);

        rotVelocity -= rotVelocity * rotFriction;
        incVelocity -= incVelocity * incFriction;

        if (RadialOrientation.Inclination > Mathf.PI / 2.1f) RadialOrientation.Inclination = Mathf.PI / 2.1f;
        if (RadialOrientation.Inclination < -Mathf.PI / 2.1f) RadialOrientation.Inclination = -Mathf.PI / 2.1f;

        if (Input.mouseScrollDelta.sqrMagnitude > 0.0f)
        {
            float scrollMagnitude = Input.mouseScrollDelta.magnitude;
            if (Input.mouseScrollDelta.y > 0.0f) // scroll up
                Distance /= 1.0f + scrollMagnitude / ScrollScale;
            if (Input.mouseScrollDelta.y < 0.0f) // scroll down
                Distance *= 1.0f + scrollMagnitude / ScrollScale;
        }
    }
	
	// Update is called once per frame
	void FixedUpdate () {
	    ProcessMouse();

	    transform.position = RadialOrientation.toCartesian.normalized * Distance;
        transform.LookAt(CenterPosition);
	}
}
