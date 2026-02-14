use crate::domain::model::NetworkBlueprint;
use cfd_core::error::Result;

pub trait GraphSink {
    type Output;

    fn build(&self, blueprint: &NetworkBlueprint) -> Result<Self::Output>;
}
