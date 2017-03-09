using System;
using UnityEngine;
using System.Linq;
using System.Collections;

public class ScaleAutoAdjust : MonoBehaviour
{

    public float Unity = 1000.0f;

    public Simulator Simulator;
    public Camera Reference;

    private float initialDrawScale;

	// Use this for initialization
	void Start ()
	{
	    Reference = FindObjectsOfType<Camera>().SingleOrDefault(c => c.enabled);
	    if (Reference == null)
	    {
	        throw new Exception("There must be exactly ONE active camera in the scene!");
	    }

	    Simulator = FindObjectsOfType<Simulator>().Single();
	    initialDrawScale = Simulator.DrawScale;
	}
	
	// Update is called once per frame
	void Update ()
	{
	    var dist = Reference.transform.position.magnitude;

        //  Make the stars larger at large distances, to increase visibility (can't SSAA 4x)
        float distScale = dist / Unity;
        distScale = Mathf.Clamp(distScale, 1.0f, 10.0f);

        Simulator.DrawScale = initialDrawScale / distScale;
    }
}
