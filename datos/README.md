# Datos

Los datos no se incluyen en el repositorio. Son públicos (licencia CC BY 4.0) y se descargan de figshare, donde están tanto los CSV como las 270 imágenes TDI:

https://figshare.com/articles/dataset/BC_cardiotox_A_cardiotoxicity_dataset_for_breast_cancer_patients/22650748

Los scripts esperan esta estructura:

```
datos/
  BC_cardiotox_clinical_variables.csv
  BC_cardiotox_functional_variable.csv
  BC_cardiotox_clinical_and_functional_variables.csv
IMAGES/
  image (1).png ... image (270).png
```

`IMAGES/` va en la raíz del repositorio, al mismo nivel que `datos/`. Las rutas se pueden cambiar en `python/config.py` y en `R/tfg_comun.R`.

Si usas los datos, cita el artículo original:

Piñeiro-Lamas, B., López-Cheda, A., Cao, R. et al. A cardiotoxicity dataset for breast cancer patients. *Sci Data* 10, 527 (2023). https://doi.org/10.1038/s41597-023-02419-1
