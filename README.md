{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "a1ded3d5-01c4-42ab-b19a-40f2bb99bdbd",
   "metadata": {},
   "source": [
    "# Subduction system response to ribbon collision: implications on the intra-plate force balance and the style of slab deformation\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "25f74fc6-a016-48d3-8b7b-01fc55ad6e76",
   "metadata": {},
   "source": [
    "[![DOI](https://zenodo.org/badge/160595955.svg)](https://zenodo.org/badge/latestdoi/160595955)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "473614de-083b-46b6-a2b3-81a8abcef84c",
   "metadata": {},
   "source": [
    "## Abstract"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "64c0029a-48cb-4306-8b72-24d349a5bc11",
   "metadata": {},
   "source": [
    "Ribbon collision is a process that can rapidly disturb the symmetry of subduction zones. Previous studies have demonstrated how ribbon collision causes rotation at the surface and contortion in the slab, but have focused on the surface kinematics. We use three-dimensional mechanical models to investigate how the dynamic evolution of ribbon collision perturbs the strain and stress field at the surface, the deformation style in the slab and the force balance in intra-plate regions. In our set of numerical simulations, we vary the angle between the trench and the ribbon to explore the slab response to ribbons colliding at different orientations. Our numerical simulations show that ribbon collision causes significant heterogeneity of stress, strain rate and vorticity in both the overriding and subducting plates surface and the slab. Slab deformation shows compartmentalization into low and high strain-rate regions around a high vorticity zone, with strain-rate variations of up to and order of magnitude occurring in both the along-strike and down-dip directions. The simulations show that changes in the collision orientation lead to contrasting styles of intra-slab deformation, where regions of high strain might yield and cause strain localization, fault activity and intra-slab seismicity. In the context of our idealised oceanic-continental subduction system, the simulations show that intra-plate stresses are affected to a similar degree by buoyancy contrasts (i.e. gravitational potential energy variations), slab-pull and ribbon collision. This partitioning allows for significant heterogeneity in the intra-plate stress regime. This work highlights how the rapid changes in strain-rate within the slab, caused by ribbon collision, can explain the seismicity gaps observed in collisional margins, which are often interpreted as slab-tears. \n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "6408080a-e33c-4c32-a419-ea9f2121b19b",
   "metadata": {},
   "source": [
    "### Surface-view"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "4a6a1d9b-06f2-4a8b-8e2d-ed6206804b74",
   "metadata": {},
   "source": [
    "<img src=\"./images/Surface_view.png\" alt=\"Drawing\" style=\"width: 1500px;/\">"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "27451d51-7c59-4c54-977f-36c09941f3f6",
   "metadata": {},
   "source": [
    "### Slab-view"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "2071684a-1d26-4fcb-b6fc-48ca749f7477",
   "metadata": {},
   "source": [
    "<img src=\"./images/Slab_view.png\" alt=\"Drawing\" style=\"width: 1500px;/\">"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "fee43cd3-0b7f-4082-8d1a-5925e076f2fe",
   "metadata": {},
   "source": [
    "## Force Balance"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "3df5e09a-1012-4fad-9714-ee21de97acae",
   "metadata": {},
   "source": [
    "<img src=\"./images/ForceBalance.png\" alt=\"Drawing\" style=\"width: 1500px;/\">"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "fa61fe92-ad9e-4b82-8d73-043dc7b5d8de",
   "metadata": {},
   "source": [
    "## This repository"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "e7ad8350-c8ae-487f-b31c-4326806dde87",
   "metadata": {},
   "source": [
    "Contains the original notebooks  scripts used to ran the numerical simulations presented in the paper: \"Subduction system response to ribbon collision: implications on the intra-plate force balance and the style of slab deformation\", within the folder \"UWGeodynamics_Scripts\". Additionally, it contains the data for computing the force balance plot. The data of the 2D models is not available due to GitHub limitations on filesizes. If need you need those files, contact andres.rodriguez1@sydney.edu.au. "
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.16"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
