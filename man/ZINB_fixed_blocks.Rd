% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ZINB_fixed_blocks-Class.R
\name{ZINB_fixed_blocks}
\alias{ZINB_fixed_blocks}
\title{R6 class for a generic normal model}
\description{
R6 class for a generic normal model

R6 class for a generic normal model
}
\section{Super class}{
\code{\link[normalblockr:NB]{normalblockr::NB}} -> \code{ZINB_fixed_blocks}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{posterior_par}}{a list with the parameters of posterior distribution W | Y}

\item{\code{entropy}}{Entropy of the conditional distribution}

\item{\code{nb_param}}{number of parameters in the model}

\item{\code{model_par}}{a list with model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix), kappa (zero-inflation probabilities)}

\item{\code{fitted}}{Y values predicted by the model}

\item{\code{who_am_I}}{a method to print what model is being fitted}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-ZINB_fixed_blocks-new}{\code{ZINB_fixed_blocks$new()}}
\item \href{#method-ZINB_fixed_blocks-clone}{\code{ZINB_fixed_blocks$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="candidates_merge"><a href='../../normalblockr/html/NB.html#method-NB-candidates_merge'><code>normalblockr::NB$candidates_merge()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="candidates_split"><a href='../../normalblockr/html/NB.html#method-NB-candidates_split'><code>normalblockr::NB$candidates_split()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="latent_network"><a href='../../normalblockr/html/NB.html#method-NB-latent_network'><code>normalblockr::NB$latent_network()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="merge"><a href='../../normalblockr/html/NB.html#method-NB-merge'><code>normalblockr::NB$merge()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="optimize"><a href='../../normalblockr/html/NB.html#method-NB-optimize'><code>normalblockr::NB$optimize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="plot_loglik"><a href='../../normalblockr/html/NB.html#method-NB-plot_loglik'><code>normalblockr::NB$plot_loglik()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="plot_network"><a href='../../normalblockr/html/NB.html#method-NB-plot_network'><code>normalblockr::NB$plot_network()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="split"><a href='../../normalblockr/html/NB.html#method-NB-split'><code>normalblockr::NB$split()</code></a></span></li>
<li><span class="pkg-link" data-pkg="normalblockr" data-topic="NB" data-id="update"><a href='../../normalblockr/html/NB.html#method-NB-update'><code>normalblockr::NB$update()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ZINB_fixed_blocks-new"></a>}}
\if{latex}{\out{\hypertarget{method-ZINB_fixed_blocks-new}{}}}
\subsection{Method \code{new()}}{
Create a new [`ZINB_fixed_blocks_fixed_sparsity`] object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ZINB_fixed_blocks$new(data, C, sparsity = 0, control = NB_control())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{object of NBData class, with responses and design matrix}

\item{\code{C}}{clustering matrix C_jq = 1 if species j belongs to cluster q}

\item{\code{sparsity}}{to apply on variance matrix when calling GLASSO}

\item{\code{control}}{structured list of more specific parameters, to generate with NB_control}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new [`ZINB_fixed_blocks_fixed_sparsity`] object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ZINB_fixed_blocks-clone"></a>}}
\if{latex}{\out{\hypertarget{method-ZINB_fixed_blocks-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ZINB_fixed_blocks$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
