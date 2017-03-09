using UnityEngine;
using System.Collections;
using UnityStandardAssets.ImageEffects;

[ExecuteInEditMode]
public class AnalyticFilters : ImageEffectBase {

	public enum Filter
	{
		None = 0,
		RedGreenOpponency,
		BlueYellowOpponency,
		WhiteBlackOpponency
	}

	public Filter FilterMode;

	void OnRenderImage(RenderTexture source, RenderTexture destination)
	{
		material.SetInt("_FilterMode", (int)FilterMode);

		Graphics.Blit(source, destination, material);
	}
}
