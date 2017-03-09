using System;
using UnityEngine;
using System.Collections;
using System.Diagnostics;
using System.Linq;

public class Simulator : MonoBehaviour
{
    private IntPtr simCtx;
    private GameObject[] objects;
    private Material[] generatedMaterials;
    private float[] sizes;

    public int Count = 500;
    public Transform BodyPrefab;
    public int StepsPerSecond = 30;
    public float DrawScale = 30.0f;

    //  https://gist.github.com/paulkaplan/5184275
    Color ColorTemperatureToRGB(float kelvin)
    {
        var temp = kelvin / 100;

        float red, green, blue;

        if (temp <= 66)
        {
            red = 255;

            green = temp;
            green = (float) (99.4708025861 * Math.Log(green) - 161.1195681661);


            if (temp <= 19)
            {
                blue = 0;
            }
            else
            {
                blue = temp - 10;
                blue = (float) (138.5177312231 * Math.Log(blue) - 305.0447927307);
            }
        }
        else
        {
            red = temp - 60;
            red = (float) ((float)329.698727446 * Math.Pow(red, -0.1332047592));

            green = temp - 60;
            green = (float) ((float)288.1221695283 * Math.Pow(green, -0.0755148492));

            blue = 255;
        }

        return new Color(
            Mathf.Clamp(red, 0.0f, 255.0f) / 255,
            Mathf.Clamp(green, 0.0f, 255.0f) / 255,
            Mathf.Clamp(blue, 0.0f, 255.0f) / 255
            );
    }

    // Use this for initialization
    void Start ()
	{
	    if (BodyPrefab == null)
	        throw new ArgumentNullException();

	    simCtx = QuickDataGravity.start_new();
        QuickDataGravity.generate_bodies(simCtx, Count);

        sizes = new float[Count];
        QuickDataGravity.get_body_radii(simCtx, sizes, Count);

	    float min = sizes.Min();
	    float max = sizes.Max();

        generatedMaterials = new Material[32];
        var refMaterial = BodyPrefab.GetComponent<Renderer>().sharedMaterial;
        for (int i = 0; i < generatedMaterials.Length; i++)
        {
            var material = new Material(refMaterial.shader);
            material.CopyPropertiesFromMaterial(refMaterial);

            float size = min + (max - min) * (i / (float)generatedMaterials.Length);
            float mass = 4.0f * Mathf.Pow(size, 2.0f);
            //float norm = (size - min) / (max - min);
            //Color c = MaxColor*norm + MinColor*(1 - norm);
            Color c = ColorTemperatureToRGB(5e8f / mass);

            material.color = c;

            generatedMaterials[i] = material;
        }
        

        objects = new GameObject[Count];
	    for (int i = 0; i < Count; i++)
	    {
	        float size = sizes[i];
            int bestMaterialIndex = (int)Math.Round((size - min) / (max - min) * (generatedMaterials.Length - 1));

            if (bestMaterialIndex >= generatedMaterials.Length)
                Debugger.Break();

            var material = generatedMaterials[bestMaterialIndex];

	        var obj = GameObject.Instantiate(BodyPrefab);
            obj.localScale = new Vector3(size, size, size) / DrawScale;
	        obj.parent = transform;
	        obj.GetComponent<Renderer>().sharedMaterial = material;
	        objects[i] = obj.gameObject;
	    }
	}
	
	// Update is called once per frame
	void FixedUpdate () {
		if (StepsPerSecond < 1)
			StepsPerSecond = 1;

		QuickDataGravity.set_sim_timestep(simCtx, 1.0f / StepsPerSecond);
        QuickDataGravity.simulate(simCtx);

	    float[] positions = new float[Count * 3]; // 3 components (x,y,z) per position
        QuickDataGravity.get_body_positions(simCtx, positions, Count);

	    for (int i = 0; i < Count; i++)
	    {
	        Vector3 newPosition = new Vector3(
                positions[i * 3 + 0],
                positions[i * 3 + 1],
                positions[i * 3 + 2]
                );

	        objects[i].transform.localScale = Vector3.one * sizes[i] / DrawScale;
	        objects[i].transform.localPosition = newPosition;
	    }
	}
}
