Shader "Custom/FullscreenDarkness" {
	Properties {
		_MainTex ("Base (RGB)", 2D) = "white" {}
	}

	SubShader {
		Pass {
			ZTest Always Cull Off ZWrite Off
					
			CGPROGRAM
			#pragma vertex vert_img
			#pragma fragment frag
			#include "UnityCG.cginc"

			uniform sampler2D _MainTex;
			uniform int _FilterMode;

			fixed4 frag (v2f_img i) : SV_Target
			{
				const fixed3 R = fixed3(1.0, 0.0, 0.0);
				const fixed3 G = fixed3(0.0, 1.0, 0.0);
				const fixed3 B = fixed3(0.0, 0.0, 1.0);
				const fixed3 Y = fixed3(1.0, 1.0, 0.0);

				fixed4 original = tex2D(_MainTex, i.uv);

				fixed3 ex = fixed3(
					(original.r - original.g) / sqrt(2.0),
					(original.r + original.g - 2.0*original.b) / sqrt(6.0),
					(original.r + original.g + original.b) / sqrt(3.0)
				);

				ex = clamp(ex, 0.0, 1.0);

				fixed4 output;

				switch (_FilterMode)
				{
				default:
					output = original;
					break;

				case (1):
					output.rgb = R * ex.x + G * (1.0 - ex.x);
					//output.rgb = R * ex.x;
					break;

				case (2):
					output.rgb = Y * ex.y + B * (1.0 - ex.y);
					break;

				case (3):
					output.rgb = fixed3(ex.z, ex.z, ex.z);
					break;
				}

				//output.rgb = ex;
				//output.rgb = fixed3(ex.x, ex.x, ex.x);
				
				output.a = original.a;
				return output;
			}
			ENDCG
		}
	}

	Fallback off
}
